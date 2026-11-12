/* ====== init ====== */
function renderAll() { renderBeam(); renderPhys(); renderAlpha(); renderRun(); renderSetup(); drawStack(); renderSummary(); renderPlots(true); }
loadLocal(); renderAll(); alibRestore(false).then(ok => { renderAlpha(); if (ok) alibAfterIndex(); });
setTimeout(() => { try { startRun(); } catch (e) { console.error(e); } }, 50);
