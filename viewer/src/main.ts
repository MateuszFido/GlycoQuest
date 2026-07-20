// Copyright (c) ETH Zurich, Mateusz Fido

import { mountViewer } from './index';
import './styles/viewer.css';

function serveHint(): string {
  if (location.protocol !== 'file:') return '';
  return (
    ' Serve this folder over HTTP, e.g. <code>cd results/viewer && python3 -m http.server 8080</code> ' +
    'then open <code>http://localhost:8080</code>.'
  );
}

const root = document.getElementById('app');
if (root && location.protocol !== 'file:') {
  root.innerHTML = '<p class="gq-empty">Loading viewer data…</p>';
  mountViewer(root, { bundleUrl: './viewer.json' }).catch((err) => {
    root.innerHTML =
      `<p class="gq-empty">Failed to load viewer: ${err}.${serveHint()}</p>`;
  });
}
