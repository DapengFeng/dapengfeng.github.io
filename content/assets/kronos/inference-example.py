import numpy as np
import pandas as pd
import torch
from model import Kronos, KronosTokenizer, KronosPredictor

cols = ["open", "high", "low", "close", "volume", "amount"]
M, lookback, horizon = 512, 400, 60
hist = pd.read_csv("history.csv", parse_dates=["timestamps"])
future = pd.read_csv("future.csv", parse_dates=["timestamps"])
assert 1 <= lookback <= M and 1 <= horizon <= M
assert len(hist) >= lookback and len(future) >= horizon
hist = hist.tail(lookback).reset_index(drop=True)
xs = hist["timestamps"]
ys = future["timestamps"].iloc[:horizon].reset_index(drop=True)
for stamps in (xs, ys):
    assert stamps.notna().all()
    assert stamps.is_monotonic_increasing and stamps.is_unique
assert ys.iloc[0] > xs.iloc[-1]
x = hist[cols].astype(float)
assert np.isfinite(x.to_numpy()).all()
assert (x["high"] >= x[["open", "close", "low"]].max(axis=1)).all()
assert (x["low"] <= x[["open", "close", "high"]].min(axis=1)).all()
assert (x[["volume", "amount"]] >= 0).all().all()

torch.manual_seed(123)
tokenizer = KronosTokenizer.from_pretrained(
    "NeoQuasar/Kronos-Tokenizer-base",
    revision="0e0117387f39004a9016484a186a908917e22426",
).eval()
model = Kronos.from_pretrained(
    "NeoQuasar/Kronos-small",
    revision="901c26c1332695a2a8f243eb2f37243a37bea320",
).eval()
predictor = KronosPredictor(model, tokenizer, device="cpu", max_context=M)
pred = predictor.predict(
    df=x, x_timestamp=xs, y_timestamp=ys, pred_len=horizon,
    T=1.0, top_k=0, top_p=0.9, sample_count=1, verbose=False,
)
assert pred.shape == (horizon, 6)
assert np.isfinite(pred.to_numpy()).all()
valid = (
    (pred["high"] >= pred[["open", "close", "low"]].max(axis=1))
    & (pred["low"] <= pred[["open", "close", "high"]].min(axis=1))
    & (pred[["volume", "amount"]] >= 0).all(axis=1)
)
print("OHLC/activity validity:", float(valid.mean()))
print(pred.head())
pred.to_csv("forecast.csv")
