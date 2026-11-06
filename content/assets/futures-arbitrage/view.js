import { calculateCalendarSpread } from './calculator.js';

const root = document.querySelector('[data-futures-calculator]');
if (root) {
  const form = root.querySelector('form');
  const error = root.querySelector('[role="alert"]');
  const format = (v, signed = false) => new Intl.NumberFormat('en-US', {
    maximumFractionDigits: 2, signDisplay: signed ? 'exceptZero' : 'auto',
  }).format(v).replace(/-/g, '−');

  function update() {
    const inputs = [...form.querySelectorAll('input[type="number"]')];
    let valid = true;
    for (const input of inputs) {
      const ok = input.value.trim() !== '' && input.checkValidity() && Number.isFinite(input.valueAsNumber);
      input.setAttribute('aria-invalid', String(!ok));
      valid &&= ok;
    }
    let result;
    if (valid) {
      try {
        result = calculateCalendarSpread({
          nearOpen: 4000, farOpen: 4100,
          ...Object.fromEntries(inputs.map(input => [input.name, input.valueAsNumber])),
          longNear: form.elements['fut-direction'].value === 'long',
        });
        valid = Object.values(result).every(Number.isFinite);
      } catch { valid = false; }
    }
    error.hidden = valid;
    for (const output of root.querySelectorAll('[data-result]')) {
      output.textContent = valid ? format(result[output.dataset.result], ['nearPnl', 'farPnl', 'net'].includes(output.dataset.result)) : '—';
    }
    const extent = valid ? Math.max(100, Math.abs(result.nearPnl), Math.abs(result.farPnl), Math.abs(result.net)) * 1.1 : 1;
    for (const bar of root.querySelectorAll('[data-bar]')) {
      bar.hidden = !valid;
      if (!valid) continue;
      const value = result[bar.dataset.bar];
      const width = Math.abs(value) / extent * 50;
      bar.style.width = `${width}%`;
      bar.style.left = `${value < 0 ? 50 - width : 50}%`;
      bar.dataset.negative = String(value < 0);
    }
  }
  for (const control of form.querySelectorAll('[disabled]')) control.disabled = false;
  form.addEventListener('input', update);
  form.addEventListener('change', update);
  form.addEventListener('submit', event => event.preventDefault());
  form.addEventListener('reset', event => {
    event.preventDefault();
    for (const input of form.querySelectorAll('input')) {
      if (input.type === 'radio') input.checked = input.defaultChecked;
      else input.value = input.defaultValue;
    }
    update();
  });
  update();
}
