#!/usr/bin/env python3
"""Recompute hypothetical futures examples with the Python standard library.

No market prices, forecasts, broker parameters or historical backtest are used.
Run: python3 calculate-examples.py
All money amounts are CNY; prices are CNY per tonne; quantity is tonnes.
"""

import json
import math


def basis_exit(spot_entry, futures_entry, spot_exit, futures_exit, quantity, cost):
    """Long physical spot and short futures, closed simultaneously."""
    spot_pnl = (spot_exit - spot_entry) * quantity
    futures_pnl = (futures_entry - futures_exit) * quantity
    entry_basis = spot_entry - futures_entry
    exit_basis = spot_exit - futures_exit
    gross = spot_pnl + futures_pnl
    assert math.isclose(gross, (exit_basis - entry_basis) * quantity)
    return {
        'spot_exit': spot_exit, 'futures_exit': futures_exit,
        'exit_basis_spot_minus_futures': exit_basis,
        'spot_pnl': spot_pnl, 'futures_pnl': futures_pnl,
        'gross_pnl': gross, 'illustrative_fixed_cost': cost, 'net_pnl': gross - cost,
    }


def residual_trade(a_entry, b_entry, a_exit, b_exit, beta, mean, std,
                   a_multiplier, b_multiplier, a_lots, b_lots, cost):
    """Short A / long B, with multipliers matching a short residual exposure."""
    exposure = a_lots * a_multiplier
    assert math.isclose(b_lots * b_multiplier, beta * exposure)
    entry_residual = a_entry - beta * b_entry
    exit_residual = a_exit - beta * b_exit
    a_pnl = -(a_exit - a_entry) * exposure
    b_pnl = (b_exit - b_entry) * b_lots * b_multiplier
    gross = a_pnl + b_pnl
    assert math.isclose(gross, exposure * (entry_residual - exit_residual))
    return {
        'a_exit': a_exit, 'b_exit': b_exit,
        'exit_residual': exit_residual, 'exit_z': (exit_residual - mean) / std,
        'short_a_pnl': a_pnl, 'long_b_pnl': b_pnl,
        'gross_pnl': gross, 'assumed_total_trading_cost': cost, 'net_pnl': gross - cost,
    }


def calculate():
    quantity = 10
    spot_entry = 2000
    futures_entry = 2060
    initial_account_cash = 3060
    margin_rate = 0.10
    settlements = [2060, 2130, 2200]
    cash_path = []
    for step, futures in enumerate(settlements):
        # Constant basis is an illustrative interim path, not expiry delivery.
        spot = futures - 60
        spot_pnl = (spot - spot_entry) * quantity
        futures_pnl = (futures_entry - futures) * quantity
        equity = initial_account_cash + futures_pnl
        required = round(futures * quantity * margin_rate)
        available = equity - required
        cash_path.append({
            'step': step, 'spot_price': spot, 'futures_settlement': futures,
            'spot_unrealized_pnl': spot_pnl,
            'futures_cumulative_pnl': futures_pnl,
            'interim_combined_price_pnl': spot_pnl + futures_pnl,
            'futures_account_equity': equity,
            'required_margin': required,
            'available_funds': available,
        })
    locked_gross_spread = (futures_entry - spot_entry) * quantity
    holding_days = 60
    annual_simple_rate = 0.0365
    spot_purchase = spot_entry * quantity
    funded_futures_and_reserve = 5000
    funded_capital = spot_purchase + funded_futures_and_reserve
    spot_interest = round(spot_purchase * annual_simple_rate * holding_days / 365, 8)
    futures_and_reserve_interest = round(funded_futures_and_reserve * annual_simple_rate * holding_days / 365, 8)
    interest = spot_interest + futures_and_reserve_interest
    storage_rate = 0.15
    storage = storage_rate * quantity * holding_days
    other_complete_costs = 60
    carrying_cost = interest + storage + other_complete_costs
    exit_scenarios = [basis_exit(spot_entry, futures_entry, s, 2100, quantity, carrying_cost)
                      for s in [2100, 2080, 2050]]
    assert [r['net_pnl'] for r in exit_scenarios] == [300, 100, -200]
    # At fixed basis, a common move in spot and futures cannot change price P&L.
    for shift in [-500, 0, 500]:
        same_basis = basis_exit(spot_entry, futures_entry, 2080 + shift, 2100 + shift,
                               quantity, carrying_cost)
        assert same_basis['net_pnl'] == 100

    initial_external_reserve = funded_futures_and_reserve - initial_account_cash
    top_up = max(0, -cash_path[-1]['available_funds'])
    account_after_top_up = cash_path[-1]['futures_account_equity'] + top_up
    remaining_external_reserve = initial_external_reserve - top_up
    delivery_price = 2100
    final_futures_pnl = (futures_entry - delivery_price) * quantity
    account_released = initial_account_cash + top_up + final_futures_pnl
    physical_proceeds = delivery_price * quantity
    closing_cash_before_costs = physical_proceeds + account_released + remaining_external_reserve
    closing_cash_after_costs = closing_cash_before_costs - carrying_cost
    final_cash_after_repayment = closing_cash_after_costs - funded_capital
    assert initial_account_cash + initial_external_reserve + spot_purchase == funded_capital
    assert (account_after_top_up, remaining_external_reserve) == (2200, 1400)
    assert (account_released, closing_cash_before_costs) == (3200, 25600)
    assert final_cash_after_repayment == locked_gross_spread - carrying_cost
    assert locked_gross_spread == 600
    assert locked_gross_spread - carrying_cost == 300
    assert [r['spot_unrealized_pnl'] for r in cash_path] == [0, 700, 1400]
    assert [r['futures_cumulative_pnl'] for r in cash_path] == [0, -700, -1400]
    assert [r['interim_combined_price_pnl'] for r in cash_path] == [0, 0, 0]
    assert [r['futures_account_equity'] for r in cash_path] == [3060, 2360, 1660]
    assert [r['required_margin'] for r in cash_path] == [2060, 2130, 2200]
    assert [r['available_funds'] for r in cash_path] == [1000, 230, -540]

    near_entry, far_entry = 4000, 4100
    near_exit, far_exit = 4050, 4110
    near_pnl = (near_exit - near_entry) * quantity
    far_pnl = (far_entry - far_exit) * quantity
    gross_pnl = near_pnl + far_pnl
    total_cost = 40
    assert (near_pnl, far_pnl, gross_pnl) == (500, -100, 400)
    assert gross_pnl == ((far_entry - near_entry) - (far_exit - near_exit)) * quantity
    assert gross_pnl - total_cost == 360

    calendar_gap_days = 30
    calendar_financing = near_entry * annual_simple_rate * calendar_gap_days / 365
    calendar_storage = 0.8 * calendar_gap_days
    calendar_other = 4
    calendar_carry = calendar_financing + calendar_storage + calendar_other
    calendar_stress_carry = calendar_carry + 20
    assert math.isclose(calendar_carry, 40)
    calendar_scenarios = [{
        'exit_far_minus_near': spread,
        'gross_pnl': ((far_entry - near_entry) - spread) * quantity,
        'net_pnl': ((far_entry - near_entry) - spread) * quantity - total_cost,
    } for spread in [40, 60, 100, 160]]
    assert [r['net_pnl'] for r in calendar_scenarios] == [560, 360, -40, -640]

    beta, residual_mean, residual_std, phi = 0.8, 60, 20, 0.9
    a_entry, b_entry = 4100, 5000
    a_multiplier, b_multiplier, a_lots, b_lots = 10, 5, 5, 8
    statistical_cost = 200
    residual_entry = a_entry - beta * b_entry
    residual_half_life = math.log(0.5) / math.log(phi)
    expected_at_half_life = residual_mean + phi ** residual_half_life * (residual_entry - residual_mean)
    statistical_scenarios = [residual_trade(
        a_entry, b_entry, a_exit, b_exit, beta, residual_mean, residual_std,
        a_multiplier, b_multiplier, a_lots, b_lots, statistical_cost,
    ) for a_exit, b_exit in [(4150, 5100), (4140, 5000)]]
    assert [r['net_pnl'] for r in statistical_scenarios] == [1300, -2200]
    assert math.isclose(expected_at_half_life, 80)
    # The fitted common component cancels only with multiplier-adjusted legs.
    for shift in [-100, 0, 100]:
        unchanged = residual_trade(
            a_entry, b_entry, a_entry + beta * shift, b_entry + shift,
            beta, residual_mean, residual_std,
            a_multiplier, b_multiplier, a_lots, b_lots, statistical_cost,
        )
        assert math.isclose(unchanged['gross_pnl'], 0)

    return {
        'kind': 'hypothetical teaching examples; not market data or a backtest',
        'units': {'money': 'CNY', 'price': 'CNY/tonne', 'quantity': 'tonnes'},
        'cash_and_carry': {
            'quantity': quantity, 'spot_entry': spot_entry,
            'futures_entry': futures_entry,
            'locked_gross_delivery_spread': locked_gross_spread,
            'assumed_complete_carrying_cost': carrying_cost,
            'assumed_delivery_net_profit': locked_gross_spread - carrying_cost,
            'initial_futures_account_cash': initial_account_cash,
            'initial_required_margin': cash_path[0]['required_margin'],
            'initial_cash_buffer': cash_path[0]['available_funds'],
            'margin_rate': margin_rate,
            'interim_path': cash_path,
            'shortfall_before_a_deposit_or_forced_close': 540,
            'basis_exits': {
                'entry_basis_spot_minus_futures': spot_entry - futures_entry,
                'scenarios': exit_scenarios,
                'notes': [
                    'Cost is held at 300 only to isolate exit-basis sensitivity; this is not an early-exit cost estimate.',
                    'An actual early exit must recalculate holding days, financing, storage, spot liquidation and trading costs.',
                    'These three marks are constructed scenarios, not forecasts or probabilities.',
                ],
            },
            'cost_breakdown': {
                'holding_days': holding_days, 'day_count_basis': 365,
                'annual_simple_financing_rate': annual_simple_rate,
                'spot_purchase_capital': spot_purchase,
                'futures_and_reserve_capital': funded_futures_and_reserve,
                'total_funded_capital': funded_capital,
                'spot_financing_cost': spot_interest,
                'futures_and_reserve_financing_cost': futures_and_reserve_interest,
                'total_financing_cost': interest,
                'storage_per_tonne_per_day': storage_rate,
                'storage_cost': storage,
                'other_complete_costs': other_complete_costs,
                'additional_net_tax_cost': 0,
                'total_cost': carrying_cost,
                'cost_per_tonne': carrying_cost / quantity,
                'delivery_break_even_futures_price': spot_entry + carrying_cost / quantity,
                'notes': [
                    'All capital is assumed borrowed for the full 60 days; idle reserve earns no interest.',
                    'Other costs are an assumed complete allowance for applicable insurance, logistics, execution and delivery.',
                    'Spot and futures use a consistent tax basis, with no additional net tax cost in this constructed example.',
                    'All expenses including financing are assumed payable at the end to make the ledger explicit.',
                ],
            },
            'funding_ledger': {
                'initial_borrowed_capital': funded_capital,
                'initial_spot_purchase': spot_purchase,
                'initial_futures_account_cash': initial_account_cash,
                'initial_external_reserve': initial_external_reserve,
                'transfer_from_reserve_to_account': top_up,
                'account_equity_after_top_up': account_after_top_up,
                'remaining_external_reserve': remaining_external_reserve,
                'assumed_delivery_price': delivery_price,
                'physical_delivery_proceeds': physical_proceeds,
                'final_futures_cumulative_pnl': final_futures_pnl,
                'futures_account_cash_released': account_released,
                'total_cash_before_costs': closing_cash_before_costs,
                'total_costs_paid_at_end': carrying_cost,
                'cash_after_costs': closing_cash_after_costs,
                'loan_principal_repaid': funded_capital,
                'net_cash_after_repayment': final_cash_after_repayment,
                'net_profit_over_capital_committed_for_60_days': final_cash_after_repayment / funded_capital,
                'notes': [
                    'The top-up is an internal cash transfer, not an additional expense or fresh borrowing.',
                    'The remaining 5000 allocation is sufficient only for this stated settlement path, not a worst-case funding ceiling.',
                    'Profit divided by total committed capital is a 60-day ratio, not ROE or a compounded annual return.',
                    'Delivery settlement accounting is summarized by physical proceeds plus cumulative futures P&L; broker details may differ.',
                ],
            },
            'notes': [
                'Interim marks hold futures minus spot at 60; they do not model expiry convergence.',
                'Interim combined price P&L of zero excludes the initially locked gross delivery spread of 600.',
                'Futures account values omit fees, interest, deposits and withdrawals.',
                'Negative available funds shows a hypothetical shortfall before a deposit or forced close.',
                'Unrealized spot gains cannot automatically fund the futures account.',
            ],
        },
        'calendar_spread': {
            'quantity_each_leg': quantity,
            'position': 'long near month, short far month',
            'entry': {'near': near_entry, 'far': far_entry, 'far_minus_near': far_entry - near_entry},
            'exit': {'near': near_exit, 'far': far_exit, 'far_minus_near': far_exit - near_exit},
            'near_leg_pnl': near_pnl, 'far_leg_pnl': far_pnl,
            'gross_pnl': gross_pnl, 'assumed_total_trading_cost': total_cost,
            'net_pnl': gross_pnl - total_cost,
            'carry_reference': {
                'delivery_gap_days': calendar_gap_days,
                'annual_simple_financing_rate': annual_simple_rate,
                'financing_per_tonne': calendar_financing,
                'storage_per_tonne_per_day': 0.8,
                'storage_per_tonne': calendar_storage,
                'other_verified_costs_per_tonne': calendar_other,
                'base_carry_per_tonne': calendar_carry,
                'stress_carry_per_tonne': calendar_stress_carry,
                'round_trip_trading_cost_per_tonne': total_cost / quantity,
                'break_even_exit_spread': far_entry - near_entry - total_cost / quantity,
                'notes': [
                    'The 40 and 60 carry values are cost scenarios, not guaranteed future spreads or a no-arbitrage interval.',
                    'A physical carry route requires eligible delivery, funding, storage and matched specifications.',
                    'Pure futures legs have no assured convergence between different delivery months.',
                    'The assumed cost reference does not include the separate 4 per tonne two-leg trading cost.',
                ],
            },
            'exit_scenarios': calendar_scenarios,
            'notes': ['Only entry and exit observations are given; connecting lines are schematic.'],
        },
        'statistical_spread': {
            'kind': 'constructed example with assumed past-window parameters; not a fitted backtest',
            'residual_definition': 'X = A - beta * B; uncentered modeled spread, not the intercept-adjusted regression residual',
            'position': 'short A, long B',
            'parameters_frozen_before_entry': {
                'beta': beta, 'mean': residual_mean, 'standard_deviation': residual_std,
                'ar1_phi': phi,
            },
            'entry': {'a': a_entry, 'b': b_entry, 'residual': residual_entry,
                      'z': (residual_entry - residual_mean) / residual_std},
            'positions': {
                'a_multiplier': a_multiplier, 'b_multiplier': b_multiplier,
                'a_short_lots': a_lots, 'b_long_lots': b_lots,
                'a_price_pnl_coefficient': -a_multiplier * a_lots,
                'b_price_pnl_coefficient': b_multiplier * b_lots,
                'short_residual_exposure': a_multiplier * a_lots,
            },
            'half_life_in_bars': residual_half_life,
            'conditional_mean_at_half_life': expected_at_half_life,
            'scenarios': statistical_scenarios,
            'break_even_exit_residual': residual_entry - statistical_cost / (a_multiplier * a_lots),
            'notes': [
                'Beta, mean, standard deviation and phi are assumed, not estimated from real prices.',
                'A and B prices and contract multipliers use compatible quote and payoff units.',
                'A fractional half-life interpolates the conditional mean; it is not a realized exit time or a guarantee.',
                'Trading costs cover both legs and both entry and exit; funding and margin stress require separate evaluation.',
                'Parameters must be estimated using only past observations, then frozen or updated under a declared out-of-sample rule.',
            ],
        },
    }


if __name__ == '__main__':
    print(json.dumps(calculate(), ensure_ascii=False, indent=2))
