#!/usr/bin/env python3
"""Recompute hypothetical futures examples with the Python standard library.

No market prices, forecasts, broker parameters or historical backtest are used.
Run: python3 calculate-examples.py
All money amounts are CNY; prices are CNY per tonne; quantity is tonnes.
"""

import json


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
    carrying_cost = 300
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
            'notes': ['Only entry and exit observations are given; connecting lines are schematic.'],
        },
    }


if __name__ == '__main__':
    print(json.dumps(calculate(), ensure_ascii=False, indent=2))
