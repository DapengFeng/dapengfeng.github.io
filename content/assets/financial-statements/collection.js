import {statements} from './model.js';
const figure = document.querySelector('[data-collection-example]');
if (figure) {
  const input = figure.querySelector('input');
  function update() {
    const values = statements(Number(input.value));
    for (const output of figure.querySelectorAll('[data-account]')) {
      output.textContent = String(values[output.dataset.account]);
    }
    figure.querySelector('[data-collected-bar]').style.width = `${values.collected}%`;
    figure.querySelector('[data-receivables-bar]').style.width = `${values.receivables}%`;
  }
  input.disabled = false;
  input.addEventListener('input', update);
  update();
}
