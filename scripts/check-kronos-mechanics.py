#!/usr/bin/env python3
"""Check Kronos mechanics offline with tiny random CPU models, not market forecasts.

Usage:
  python scripts/check-kronos-mechanics.py --repo /path/to/reviewed/Kronos \
    --output content/assets/kronos/mechanics-results.json

Use CPU PyTorch 2.10.0 and numpy==2.2.6, pandas==2.2.2, einops==0.8.1,
huggingface_hub==0.33.1, safetensors==0.6.2, tqdm==4.67.1.
Dependencies and the pinned upstream checkout must already exist. This script
never installs packages or downloads weights. Temporary random checkpoints are
removed automatically. Omit --output to write /tmp/kronos-mechanics-results.json.
"""
import argparse
import json
import os
import platform
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory

SHA = '67b630e67f6a18c9e9be918d9b4337c960db1e9a'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', required=True, type=Path,
                        help='Reviewed Kronos Git checkout at the pinned source SHA.')
    parser.add_argument('--output', type=Path, default=Path('/tmp/kronos-mechanics-results.json'),
                        help='JSON result path (default: %(default)s).')
    args = parser.parse_args()
    if not __debug__:
        parser.error('Assertions must be enabled; run Python without -O or -OO.')
    repo = args.repo.expanduser().resolve()
    try:
        actual_sha = subprocess.check_output(
            ['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True, stderr=subprocess.PIPE
        ).strip()
    except (OSError, subprocess.CalledProcessError) as exc:
        parser.error(f'Cannot read the upstream Git checkout: {exc}')
    if actual_sha != SHA:
        parser.error(f'Expected upstream commit {SHA}; found {actual_sha}.')
    if subprocess.run(['git', '-C', str(repo), 'diff', '--quiet', 'HEAD', '--', 'model'],
                      check=False).returncode:
        parser.error('Upstream model source has local changes; use the reviewed revision.')

    os.environ['HF_HUB_OFFLINE'] = '1'
    os.environ['HF_HUB_DISABLE_TELEMETRY'] = '1'
    sys.path.insert(0, str(repo))
    import numpy as np
    import pandas as pd
    import torch
    from model import Kronos, KronosTokenizer, KronosPredictor
    from model.module import MultiHeadCrossAttentionWithRoPE
    from model.kronos import auto_regressive_inference, calc_time_stamps

    torch.manual_seed(20261026)
    np.random.seed(20261026)
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)

    result = {
        'scope': 'Offline mechanics only; tiny randomly initialized CPU models; no pretrained weights or financial accuracy measured.',
        'source_sha': SHA,
        'seed': 20261026,
        'environment': {
            'python': platform.python_version(), 'platform': platform.platform(),
            'torch': torch.__version__, 'device': 'cpu', 'threads': torch.get_num_threads(),
            'packages': {p: version(p) for p in ['numpy', 'pandas', 'einops', 'huggingface_hub', 'safetensors', 'tqdm']},
        },
        'checks': {},
    }
    C = result['checks']

    def diff(a, b):
        return float((a-b).abs().max())

    tok = KronosTokenizer(
        d_in=6, d_model=16, n_heads=4, ff_dim=32,
        n_enc_layers=4, n_dec_layers=4,
        ffn_dropout_p=0, attn_dropout_p=0, resid_dropout_p=0,
        s1_bits=2, s2_bits=2, beta=.05, gamma0=1., gamma=1.1, zeta=.05, group_size=2,
    ).eval()
    model = Kronos(
        s1_bits=2, s2_bits=2, n_layers=2, d_model=16, n_heads=4, ff_dim=32,
        ffn_dropout_p=0, attn_dropout_p=0, resid_dropout_p=0, token_dropout_p=0, learn_te=True,
    ).eval()
    result['tiny_config'] = {'B': 2, 'L': 6, 'features': 6, 'bits': [2, 2], 'd_model': 16, 'n_heads': 4,
                             'tokenizer_config_layers': [4, 4], 'tokenizer_actual_blocks': [len(tok.encoder), len(tok.decoder)],
                             'backbone_blocks': 2, 'dropout': 0}

    with torch.no_grad():
        x = torch.randn(2, 6, 6)
        x_changed = x.clone()
        x_changed[:, 4:] += 30*torch.randn_like(x_changed[:, 4:])
        (prefix_recon, full_recon), _, quantized, ids_full = tok(x)
        pair = tok.encode(x, half=True)
        pair_changed = tok.encode(x_changed, half=True)
        changed_recons, _, _, _ = tok(x_changed)
        assert all(torch.equal(a[:, :4], b[:, :4]) for a, b in zip(pair, pair_changed))
        assert diff(prefix_recon[:, :4], changed_recons[0][:, :4]) < 1e-6
        assert diff(full_recon[:, :4], changed_recons[1][:, :4]) < 1e-6
        C['tokenizer_causal_prefix'] = {
            'prefix_length': 4, 'suffix_changed': 2, 'prefix_ids_exactly_equal': True,
            'prefix_reconstruction_max_abs_diff': diff(prefix_recon[:, :4], changed_recons[0][:, :4]),
            'full_reconstruction_max_abs_diff': diff(full_recon[:, :4], changed_recons[1][:, :4]),
        }

        unpacked = tok.indices_to_bits(pair, half=True)
        assert torch.equal(unpacked, quantized)
        assert torch.equal(ids_full, pair[0] + 4*pair[1])
        decoded = tok.decode(pair, half=True)
        assert diff(decoded, full_recon) < 1e-6
        allpairs = torch.cartesian_prod(torch.arange(4), torch.arange(4))
        exhaustive_pair = [allpairs[:, 0][None, :], allpairs[:, 1][None, :]]
        allbits = tok.indices_to_bits(exhaustive_pair, half=True)
        assert torch.equal(tok.tokenizer.bits_to_indices(allbits[:, :, :2]), exhaustive_pair[0])
        assert torch.equal(tok.tokenizer.bits_to_indices(allbits[:, :, 2:]), exhaustive_pair[1])
        composite = tok.tokenizer.bits_to_indices(allbits)
        split_by_embedding = model.embedding.split_token(composite, 2)
        assert torch.equal(split_by_embedding[0], exhaustive_pair[1])
        assert torch.equal(split_by_embedding[1], exhaustive_pair[0])
        C['bit_roundtrip'] = {
            'all_16_pairs_roundtrip_exact': True,
            'quantized_equals_decode_unpacked_bits': True,
            'full_id_relation': 's1 + 4*s2',
            'decode_vs_forward_max_abs_diff': diff(decoded, full_recon),
            'optional_composite_embedding_route_swaps_halves': True,
            'quantized_unique_values': torch.unique(quantized).tolist(),
        }

        stamp = torch.tensor([[[i, 10, 1, 20, 3] for i in range(6)]]).repeat(2, 1, 1)
        s1_logits, context = model.decode_s1(pair[0], pair[1], stamp)
        s2_logits = model.decode_s2(context, pair[0])
        single_s2_logits = model.decode_s2(context, pair[0][:, -1:])
        logits = model(pair[0], pair[1], stamp)
        pair_mod = [t.clone() for t in pair]
        pair_mod[0][:, 4:] = (pair_mod[0][:, 4:] + 1) % 4
        pair_mod[1][:, 4:] = (pair_mod[1][:, 4:] + 2) % 4
        stamp_mod = stamp.clone()
        stamp_mod[:, 4:, 0] += 20
        other_s1, other_context = model.decode_s1(pair_mod[0], pair_mod[1], stamp_mod)
        assert diff(context[:, :4], other_context[:, :4]) < 1e-6
        changed_decoded = tok.decode(pair_mod, half=True)
        assert diff(decoded[:, :4], changed_decoded[:, :4]) < 1e-6
        C['shapes_and_backbone_causality'] = {
            'input': list(x.shape), 'quantized': list(quantized.shape), 'pair_ids': [list(t.shape) for t in pair],
            'tokenizer_reconstructions': [list(prefix_recon.shape), list(full_recon.shape)],
            'forward_logits': [list(t.shape) for t in logits],
            'context': list(context.shape), 'decode_s2_single_query_result': list(single_s2_logits.shape),
            'backbone_prefix_max_abs_diff': diff(context[:, :4], other_context[:, :4]),
            's1_prefix_logits_max_abs_diff': diff(s1_logits[:, :4], other_s1[:, :4]),
            'decoder_prefix_max_abs_diff': diff(decoded[:, :4], changed_decoded[:, :4]),
        }

        z = torch.nn.functional.normalize(torch.randn(2, 6, 4), dim=-1)
        h_sample, h_marginal, avg_prob = tok.tokenizer.bsq.soft_entropy_loss(z)
        binary_p = torch.sigmoid(-4*z/2)
        expected_avg = torch.stack([binary_p, 1-binary_p], -1).mean((0, 1))
        expected_entropy = -(expected_avg * torch.log(expected_avg+1e-8)).sum()
        assert avg_prob.shape == (4, 2)
        assert diff(avg_prob, expected_avg) < 1e-6
        C['analytical_entropy_branch'] = {
            'actual_avg_prob_shape': list(avg_prob.shape), 'group_code_probability_shape_would_be': [2, 4],
            'avg_prob_vs_per_bit_reference_max_abs_diff': diff(avg_prob, expected_avg),
            'marginal_entropy_vs_per_bit_reference_abs_diff': float(abs(h_marginal-expected_entropy)),
            'sample_entropy': float(h_sample), 'marginal_entropy': float(h_marginal),
        }

        attention = MultiHeadCrossAttentionWithRoPE(16, 4, 0., 0.)
        q = torch.randn(2, 1, 16)
        h = torch.randn(2, 6, 16)
        h_mod = h.clone()
        h_mod[:, 1:] += torch.randn_like(h_mod[:, 1:])*7
        attention.train()
        train_out = attention(q, h, h)
        train_out_changed = attention(q, h_mod, h_mod)
        first_key_reference = attention.out_proj(attention.v_proj(h[:, :1]))
        assert diff(train_out, first_key_reference) < 1e-6
        assert diff(train_out, train_out_changed) == 0
        attention.eval()
        eval_out = attention(q, h, h)
        eval_out_changed = attention(q, h_mod, h_mod)
        assert diff(eval_out, eval_out_changed) > 1e-3
        C['single_query_train_vs_eval_attention'] = {
            'dropout': 0, 'query_length': 1, 'key_length': 6,
            'train_vs_first_key_only_reference_max_abs_diff': diff(train_out, first_key_reference),
            'train_change_when_keys_1_to_5_perturbed': diff(train_out, train_out_changed),
            'eval_change_when_keys_1_to_5_perturbed': diff(eval_out, eval_out_changed),
            'train_vs_eval_max_abs_diff': diff(train_out, eval_out),
            'inference_rope_cos': attention.rotary.cos_cached.flatten().tolist(),
            'inference_rope_sin': attention.rotary.sin_cached.flatten().tolist(),
        }

        # Fix sibling queries to isolate the mask from random first-token draws.
        model.eval()
        eval_s2 = model.decode_s2(context, pair[0])
        eval_s2_changed = model.decode_s2(other_context, pair[0])
        model.train()
        train_s2 = model.decode_s2(context, pair[0])
        train_s2_changed = model.decode_s2(other_context, pair[0])
        assert diff(train_s2[:, :4], train_s2_changed[:, :4]) < 1e-6
        assert diff(eval_s2[:, :4], eval_s2_changed[:, :4]) > 1e-4
        model.eval()
        C['full_sequence_s2_prefix_dependence'] = {
            'fixed_sibling_queries': True,
            'train_prefix_change_from_future_context': diff(train_s2[:, :4], train_s2_changed[:, :4]),
            'eval_prefix_change_from_future_context': diff(eval_s2[:, :4], eval_s2_changed[:, :4]),
        }

        hist_ts = pd.Series(pd.date_range('2026-01-01 00:00', periods=6, freq='5min'))
        future_ts = pd.Series(pd.date_range('2026-01-01 00:30', periods=9, freq='5min'))
        frame = pd.DataFrame(x[0].numpy(), columns=['open', 'high', 'low', 'close', 'volume', 'amount'])
        predictor = KronosPredictor(model, tok, device='cpu', max_context=8)
        seen_stamp_last = []
        hook = model.time_emb.register_forward_pre_hook(lambda m, args: seen_stamp_last.append(args[0][:, -1, :].clone()))
        short_prediction = predictor.predict(frame, hist_ts, future_ts[:3], pred_len=3, top_k=1, top_p=1., sample_count=2, verbose=False)
        hook.remove()
        assert short_prediction.shape == (3, 6)
        assert [int(s[0, 0]) for s in seen_stamp_last] == [25, 30, 35]
        cap_output = auto_regressive_inference(
            tok, model, x[:1], torch.from_numpy(calc_time_stamps(hist_ts).values)[None],
            torch.from_numpy(calc_time_stamps(future_ts).values)[None],
            max_context=8, pred_len=9, top_k=1, top_p=1., sample_count=2,
        )
        assert cap_output.shape == (1, 8, 6)
        horizon_error = None
        try:
            predictor.predict(frame, hist_ts, future_ts, pred_len=9, top_k=1, top_p=1., verbose=False)
        except ValueError as exc:
            horizon_error = str(exc)
        assert horizon_error is not None
        C['generation_horizon_and_timestamp_alignment'] = {
            'history_length': 6, 'max_context': 8, 'normal_forecast_shape': list(short_prediction.shape),
            'last_consumed_minutes_each_step': [int(s[0, 0]) for s in seen_stamp_last],
            'target_minutes': future_ts[:3].dt.minute.tolist(),
            'pred_len_9_low_level_shape': list(cap_output.shape),
            'pred_len_9_predict_error': horizon_error,
        }

    # Test the pinned loader with our own tiny random checkpoint; never contact Hub.
    with TemporaryDirectory(prefix='kronos-mechanics-') as checkpoint_dir:
        model.eval()
        model.save_pretrained(checkpoint_dir)
        loaded = Kronos.from_pretrained(checkpoint_dir, local_files_only=True)
        assert loaded.training is True
        loaded_predictor = KronosPredictor(loaded, tok, device='cpu', max_context=8)
        assert loaded.training is True
        loaded.eval()
        assert loaded.training is False
    C['local_safetensors_eval_mode'] = {
        'saved_checkpoint_origin': 'locally initialized tiny random model',
        'training_after_from_pretrained': True,
        'training_after_predictor_constructor': True,
        'training_after_explicit_eval': False,
    }

    result['all_assertions_passed'] = True
    out = args.output.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
