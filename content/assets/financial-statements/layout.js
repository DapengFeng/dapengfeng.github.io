const article = document.querySelector('.fs-article');
const contents = article?.closest('.article-page')?.querySelector('.article-toc details');

if (contents) {
  const compact = matchMedia('(max-width: 1100px)');
  const adapt = () => { contents.open = !compact.matches; };
  adapt();
  compact.addEventListener('change', adapt);
}
