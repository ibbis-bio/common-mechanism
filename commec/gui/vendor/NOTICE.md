# Third-party notices

This directory contains a vendored third-party asset, bundled so the commec
screening report renders offline (no CDN).

## plotly.js

- File: `plotly-3.0.1.min.js`
- Version: 3.0.1
- Copyright: 2012-2025, Plotly, Inc.
- License: **MIT**
- Upstream: https://github.com/plotly/plotly.js

The MIT copyright + permission notice is retained in the file's own header
comment, and the file is shipped **unmodified** (the GUI only rewrites the
`<script src=...>` reference in commec's generated HTML to point at this local
copy: it never edits this file), so MIT's attribution requirement is met by
distributing the file as-is.

commec-gui itself is MIT-licensed (see the repository `LICENSE`, (c) IBBIS),
so bundling this MIT dependency is license-compatible. When packaging or
upstreaming, include this notice in the build's third-party / NOTICES listing.
