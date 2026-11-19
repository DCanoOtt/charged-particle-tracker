/* ====== Saving a table of the interface as PNG, PDF, Excel (.xlsx) or CSV ====== */
// tblBar(id, name) = row of buttons placed right after the <div class="tscroll"> that holds the table
const tblBar = (id, name) =>
  '<div class="btns" data-tx="' +
  id +
  '" data-txn="' +
  name +
  '" style="margin:4px 0 8px;align-items:center"><span class="hint" style="margin:0">Save this table:</span>' +
  ['png', 'pdf', 'xlsx', 'csv'].map(k => '<button class="btn" data-txk="' + k + '">' + (k === 'xlsx' ? 'Excel' : k.toUpperCase()) + '</button>').join('') +
  '<span class="hint" id="txm_' +
  id +
  '"></span></div>';
// text of the cells of a table, row by row
function tblRows(t) {
  return [...t.rows].map(r => [...r.cells].map(c => c.textContent.replace(/\s+/g, ' ').trim()));
}
const tblNum = s => (/^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(s) ? +s : null);
// "± 12" -> 12 (the column header already says ±)
const tblVal = s => s.replace(/^±\s*/, '');
function tblCsv(rows) {
  const q = s => (/[",;\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s);
  return '﻿' + rows.map((r, i) => r.map(c => q(i ? tblVal(c) : c)).join(',')).join('\r\n') + '\r\n';
}
const tblXml = s =>
  s
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/[\u0000-\u0008\u000b\u000c\u000e-\u001f]/g, '');
function tblColName(i) {
  let s = '';
  for (i++; i > 0; i = Math.floor((i - 1) / 26)) s = String.fromCharCode(65 + ((i - 1) % 26)) + s;
  return s;
}
// minimal .xlsx: numbers stay numbers, everything else is text; header row in bold
function tblXlsx(rows, sheet) {
  const enc = new TextEncoder(),
    f = (name, txt) => ({ name, data: enc.encode(txt) }),
    head = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n',
    ns = 'http://schemas.openxmlformats.org/',
    sd = rows
      .map(
        (r, i) =>
          '<row r="' +
          (i + 1) +
          '">' +
          r
            .map((c, j) => {
              const ref = tblColName(j) + (i + 1),
                v = i ? tblVal(c) : c,
                n = i ? tblNum(v) : null;
              if (n !== null && Number.isFinite(n)) return '<c r="' + ref + '"><v>' + n + '</v></c>';
              return '<c r="' + ref + '" t="inlineStr"' + (i ? '' : ' s="1"') + '><is><t xml:space="preserve">' + tblXml(v) + '</t></is></c>';
            })
            .join('') +
          '</row>'
      )
      .join(''),
    cols = (rows[0] || [])
      .map((c, j) => '<col min="' + (j + 1) + '" max="' + (j + 1) + '" width="' + Math.min(40, Math.max(12, c.length * 0.9)) + '" customWidth="1"/>')
      .join('');
  return zipMulti([
    f(
      '[Content_Types].xml',
      head +
        '<Types xmlns="' +
        ns +
        'package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/></Types>'
    ),
    f(
      '_rels/.rels',
      head +
        '<Relationships xmlns="' +
        ns +
        'package/2006/relationships"><Relationship Id="rId1" Type="' +
        ns +
        'officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'
    ),
    f(
      'xl/workbook.xml',
      head +
        '<workbook xmlns="' +
        ns +
        'spreadsheetml/2006/main" xmlns:r="' +
        ns +
        'officeDocument/2006/relationships"><sheets><sheet name="' +
        tblXml(sheet.slice(0, 31).replace(/[\\/?*[\]:]/g, '_')) +
        '" sheetId="1" r:id="rId1"/></sheets></workbook>'
    ),
    f(
      'xl/_rels/workbook.xml.rels',
      head +
        '<Relationships xmlns="' +
        ns +
        'package/2006/relationships"><Relationship Id="rId1" Type="' +
        ns +
        'officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="' +
        ns +
        'officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>'
    ),
    f(
      'xl/styles.xml',
      head +
        '<styleSheet xmlns="' +
        ns +
        'spreadsheetml/2006/main"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>'
    ),
    f(
      'xl/worksheets/sheet1.xml',
      head + '<worksheet xmlns="' + ns + 'spreadsheetml/2006/main"><cols>' + cols + '</cols><sheetData>' + sd + '</sheetData></worksheet>'
    )
  ]);
}
// table drawn on a canvas (white background, black text) for PNG and PDF
function tblCanvas(rows, title) {
  const S = 2,
    fs = 13,
    pad = 6,
    lh = 16,
    maxW = 150,
    cv = document.createElement('canvas'),
    cx = cv.getContext('2d'),
    font = b => (b ? 'bold ' : '') + fs + 'px Inter, Arial, sans-serif';
  cx.font = font(false);
  const wrap = (s, w, bold) => {
    cx.font = font(bold);
    const out = [];
    let cur = '';
    for (const word of s.split(' ')) {
      const t = cur ? cur + ' ' + word : word;
      if (cur && cx.measureText(t).width > w) {
        out.push(cur);
        cur = word;
      } else cur = t;
    }
    out.push(cur);
    return out;
  };
  const nc = Math.max(...rows.map(r => r.length)),
    cw = [];
  for (let j = 0; j < nc; j++) {
    let w = 40;
    rows.forEach((r, i) => {
      cx.font = font(i === 0);
      const t = r[j] || '';
      w = Math.max(w, Math.min(i === 0 ? maxW : 1e9, cx.measureText(t).width));
    });
    cw.push(Math.ceil(w) + 2 * pad);
  }
  const cells = rows.map((r, i) => r.map((c, j) => (i === 0 ? wrap(c, cw[j] - 2 * pad, true) : [c]))),
    rh = cells.map(r => Math.max(...r.map(l => l.length)) * lh + 2 * pad),
    th = title ? 26 : 0,
    W = cw.reduce((a, b) => a + b, 0) + 1,
    H = rh.reduce((a, b) => a + b, 0) + th + 1;
  cv.width = Math.ceil(W * S);
  cv.height = Math.ceil(H * S);
  cx.scale(S, S);
  cx.fillStyle = '#fff';
  cx.fillRect(0, 0, W, H);
  cx.textBaseline = 'top';
  if (title) {
    cx.fillStyle = '#000';
    cx.font = font(true);
    cx.fillText(title, 2, 5);
  }
  let y = th;
  cells.forEach((r, i) => {
    let x = 0;
    if (i === 0) {
      cx.fillStyle = '#e8ebf0';
      cx.fillRect(0, y, W, rh[i]);
    }
    r.forEach((lines, j) => {
      cx.fillStyle = '#000';
      cx.font = font(i === 0);
      lines.forEach((l, k) => cx.fillText(l, x + pad, y + pad + k * lh));
      cx.strokeStyle = '#9aa3b2';
      cx.lineWidth = 1;
      cx.strokeRect(x + 0.5, y + 0.5, cw[j], rh[i]);
      x += cw[j];
    });
    y += rh[i];
  });
  return cv;
}
document.addEventListener('click', async ev => {
  const b = ev.target.closest && ev.target.closest('[data-txk]');
  if (!b) return;
  const bar = b.closest('[data-tx]'),
    box = bar.previousElementSibling,
    t = box && box.querySelector('table');
  if (!t) return;
  const kind = b.dataset.txk,
    sel = '#txm_' + bar.dataset.tx,
    m = $(sel),
    rows = tblRows(t),
    name = 'cpt_' + bar.dataset.txn;
  m.textContent = 'Preparing…';
  try {
    if (kind === 'csv') await saveFile(name + '.csv', new Blob([tblCsv(rows)], { type: 'text/csv;charset=utf-8' }), 'CSV', sel);
    else if (kind === 'xlsx')
      await saveFile(
        name + '.xlsx',
        new Blob([tblXlsx(rows, bar.dataset.txn.replace(/_/g, ' '))], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' }),
        'Excel',
        sel
      );
    else {
      const cv = tblCanvas(rows, bar.dataset.txn.replace(/_/g, ' '));
      if (kind === 'png') await saveFile(name + '.png', await canvasBlob(cv, 'image/png'), 'PNG', sel);
      else await saveFile(name + '.pdf', await makePdf(cv, null), 'PDF', sel);
    }
  } catch (e) {
    m.textContent = 'Could not save (' + ((e && (e.code || e.message)) || e) + ').';
  }
});
