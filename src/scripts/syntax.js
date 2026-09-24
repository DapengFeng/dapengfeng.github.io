/* Small, dependency-free lexer. Highlighted HTML is presentation only; textarea owns the source. */
(() => {
  const words = text => new Set(text.split(' '));
  const keywords = {
    'c++': words('alignas alignof asm auto break case catch class concept const consteval constexpr constinit const_cast continue co_await co_return co_yield decltype default delete do dynamic_cast else enum explicit export extern false for friend goto if inline mutable namespace new noexcept nullptr operator private protected public register reinterpret_cast requires return sizeof static static_assert static_cast struct switch template this thread_local throw true try typedef typeid typename union using virtual volatile while'),
    rust: words('as async await break const continue crate dyn else enum extern false fn for if impl in let loop match mod move mut pub ref return self Self static struct super trait true type unsafe use where while'),
    python: words('and as assert async await break class continue def del elif else except False finally for from global if import in is lambda None nonlocal not or pass raise return True try while with yield')
  };
  const types = words('bool char char8_t char16_t char32_t double float int long short signed unsigned void wchar_t size_t ptrdiff_t string string_view vector array complex tuple pair optional span auto u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128 isize f32 f64 str String Vec Option Result Some None Ok Err');
  const escape = text => text.replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function highlight(source, language) {
    language = ({cpp:'c++',c:'c++',py:'python',rs:'rust'})[language] || language;
    if (!keywords[language]) return escape(source);
    const python = language === 'python';
    // Sticky expressions consume exactly the current position; source is never rewritten.
    const token = python
      ? /#[^\n]*|(?:[rubf]{1,2})?(?:"""[\s\S]*?(?:"""|$)|'''[\s\S]*?(?:'''|$)|"(?:\\[\s\S]|[^"\\\n])*(?:"|(?=\n|$))|'(?:\\[\s\S]|[^'\\\n])*(?:'|(?=\n|$)))|(?:0[xob][\da-f_]+|\d[\d_]*(?:\.\d[\d_]*)?(?:e[+-]?[\d_]+)?j?)\b|[A-Za-z_$][\w$]*|\s+|./giy
      : /\/\/[^\n]*|"(?:\\[\s\S]|[^"\\\n])*(?:"|(?=\n|$))|'(?:\\(?:x[\da-f]{2}|u\{[\da-f]+\}|u[\da-f]{4}|.)|[^'\\\n])'|(?:0[xb][\da-f_']+|\d[\d_']*(?:\.\d[\d_']*)?(?:e[+-]?[\d_']+)?)(?:[uifl](?:8|16|32|64|128|size)?|ull|ll|ul|lu|l|f)?\b|#[ \t]*[A-Za-z_]\w*|[A-Za-z_$][\w$]*|\s+|./giy;
    const parts = [];
    let position = 0;
    while (position < source.length) {
      let value, kind = '';
      if (!python && source.startsWith('/*', position)) {
        let end = position + 2, depth = 1;
        while (end < source.length && depth) {
          if (language === 'rust' && source.startsWith('/*', end)) { depth++; end += 2; }
          else if (source.startsWith('*/', end)) { depth--; end += 2; }
          else end++;
        }
        value = source.slice(position, end); kind = 'comment';
      } else {
        const rest = source.slice(position);
        const raw = language === 'rust' ? /^(?:br|r)(#{0,255})"/.exec(rest) : language === 'c++' ? /^(?:u8|u|U|L)?R"([^\s()\\]{0,16})\(/.exec(rest) : null;
        if (raw) {
          const close = language === 'rust' ? '"' + raw[1] : ')' + raw[1] + '"';
          const end = source.indexOf(close, position + raw[0].length);
          value = source.slice(position, end < 0 ? source.length : end + close.length); kind = 'string';
        } else {
          token.lastIndex = position;
          value = token.exec(source)?.[0] || source[position];
          if (value.startsWith('//') || (python && value.startsWith('#'))) kind = 'comment';
          else if (/^(?:[rubf]{0,2})?["']/i.test(value) && value.length > 1) kind = 'string';
          else if (/^\d/.test(value)) kind = 'number';
          else if (value.startsWith('#')) kind = 'directive';
          else if (keywords[language].has(value)) kind = 'keyword';
          else if (types.has(value)) kind = 'type';
          else if (/^[A-Za-z_$]/.test(value) && /^\s*!?\s*\(/.test(source.slice(position + value.length))) kind = 'function';
          else if (/^[A-Z][A-Z_\d]+$/.test(value)) kind = 'constant';
        }
      }
      parts.push(kind ? `<span class="syntax-${kind}">${escape(value)}</span>` : escape(value));
      position += value.length;
    }
    return parts.join('');
  }
  function mount(editor, language) {
    const surface = document.createElement('div'); surface.className = 'compiler-editor-surface';
    const gutter = document.createElement('div'); gutter.className = 'compiler-gutter'; gutter.setAttribute('aria-hidden', 'true');
    const numbers = document.createElement('div'); numbers.className = 'compiler-line-numbers'; gutter.append(numbers);
    const active = document.createElement('div'); active.className = 'compiler-active-line'; active.setAttribute('aria-hidden', 'true');
    const paint = document.createElement('pre'); paint.className = 'compiler-highlight'; paint.setAttribute('aria-hidden', 'true');
    const code = document.createElement('code'); paint.append(code);
    editor.before(surface); surface.append(active, paint, gutter, editor);
    let activeNumber;
    function selection() {
      const end = editor.selectionDirection === 'backward' ? editor.selectionStart : editor.selectionEnd;
      const line = editor.value.slice(0, end).split('\n').length - 1;
      const style = getComputedStyle(editor);
      active.style.top = `${editor.clientTop + parseFloat(style.paddingTop) + line * parseFloat(style.lineHeight) - editor.scrollTop}px`;
      activeNumber?.classList.remove('is-current');
      activeNumber = numbers.children[line]; activeNumber?.classList.add('is-current');
    }
    function scroll() {
      paint.scrollTop = editor.scrollTop; paint.scrollLeft = editor.scrollLeft;
      numbers.style.transform = `translateY(${-editor.scrollTop}px)`; selection();
    }
    function size() { paint.style.width = `${editor.clientWidth}px`; paint.style.height = `${editor.clientHeight}px`; scroll(); }
    function refresh() {
      code.innerHTML = highlight(editor.value, language) + '\n';
      const count = editor.value.split('\n').length;
      surface.style.setProperty('--editor-gutter', `${Math.max(2, String(count).length)}ch`);
      numbers.innerHTML = Array.from({length:count}, (_, i) => `<span>${i + 1}</span>`).join('');
      scroll();
    }
    editor.addEventListener('scroll', scroll, {passive:true});
    for (const event of ['select', 'click', 'keyup', 'focus']) editor.addEventListener(event, selection);
    document.addEventListener('selectionchange', () => { if (document.activeElement === editor) selection(); });
    new ResizeObserver(size).observe(editor);
    refresh(); size();
    return refresh;
  }
  window.CodeSyntax = {highlight, mount};
  document.querySelectorAll('pre > code[class*="language-"]:not([data-godbolt])').forEach(code => {
    const language = [...code.classList].find(c => c.startsWith('language-'))?.slice(9);
    code.innerHTML = highlight(code.textContent, language);
    code.parentElement.classList.add('syntax-block');
  });
})();
