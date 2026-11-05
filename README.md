# Symmetric-state research catalogues

[Open the catalogues](https://aenictusGitHub.github.io/symmetric-state-catalogues-public/).

- [Sector lengths](https://aenictusGitHub.github.io/symmetric-state-catalogues-public/symmetric-sector-catalogue/)
- [Anticoherent states](https://aenictusGitHub.github.io/symmetric-state-catalogues-public/anticoherent-state-catalogue/)
- [Anticoherent subspaces](https://aenictusGitHub.github.io/symmetric-state-catalogues-public/ac-subspace-catalogue/)
- [Majorana degeneracy](https://aenictusGitHub.github.io/symmetric-state-catalogues-public/majorana-degeneracy-catalogue/)

The catalogues distinguish proved results, exact attained candidates, numerical examples and unresolved bounds. Consult each entry's sources and certificates for the scope of its claims.

## Viewing locally

Clone or download this repository, then open `index.html` in a browser. Alternatively run `python3 -m http.server 8000` from the repository directory and visit `http://localhost:8000`.

## Publication

GitHub Actions builds and deploys the four catalogues to GitHub Pages. Full research files and exact certificates remain in this repository. To keep the published site within GitHub Pages' size limit, certificate files larger than 10 MiB are served from the corresponding commit on GitHub rather than duplicated in the site. The build rewrites their download links automatically.

Run `python3 scripts/build_site.py` to generate `_site`. No scientific data is changed by the build.
