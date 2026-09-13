/* 预渲染 KaTeX：读取 HTML 中的 \x01Mn\x02 占位符，替换为 KaTeX 的静态 HTML。
 * 用法: node render_math.js <in.html> <math.json> <out.html> <katex.min.js 路径>
 */
const fs = require('fs');

const [inHtml, mathJson, outHtml, katexPath] = process.argv.slice(2);
const katex = require(katexPath);

const html = fs.readFileSync(inHtml, 'utf8');
const items = JSON.parse(fs.readFileSync(mathJson, 'utf8'));

let nOk = 0;
const errors = [];
const out = html.replace(/\x01M(\d+)\x02/g, (m, idx) => {
  const it = items[Number(idx)];
  if (!it) return m;
  try {
    const r = katex.renderToString(it.tex, {
      displayMode: !!it.display,
      throwOnError: true,
      strict: false,
      trust: true,
    });
    nOk++;
    return r;
  } catch (e) {
    errors.push('#' + idx + ' ' + it.tex + '  ==> ' + e.message);
    return '<span style="color:#c00;border:1px solid #c00;padding:0 2px">' +
           it.tex.replace(/[<>&]/g, '') + '</span>';
  }
});

fs.writeFileSync(outHtml, out, 'utf8');
console.log('KaTeX rendered:', nOk, '/', items.length);
if (errors.length) {
  console.log('--- ERRORS (' + errors.length + ') ---');
  errors.forEach(e => console.log(e));
  process.exitCode = 2;
}
