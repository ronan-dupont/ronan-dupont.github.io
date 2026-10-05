# Ronan Dupont — Academic Website

Official personal academic website of **Ronan Dupont**, postdoctoral researcher at **Nagoya University** working on numerical linear algebra, iterative methods, sparse linear systems and scientific computing.

**Live website:** https://ronan-dupont.github.io/

Research topics include numerical linear algebra, Krylov subspace methods, preconditioning, Virtual Element Methods, Helmholtz problems and coastal morphodynamics.

Academic identifiers and profiles:

- ORCID: https://orcid.org/0000-0003-2260-0755
- J-GLOBAL: https://jglobal.jst.go.jp/en/detail?JGLOBAL_ID=202601003611729927
- ResearchGate: https://www.researchgate.net/profile/Ronan-Dupont
- Nagoya University / Zhang–Sogabe Laboratory: https://na.nuap.nagoya-u.ac.jp/member/

The site includes:

- English, French and Japanese landing pages.
- Publication pages with DOI and Google Scholar-compatible citation metadata.
- Research talks, teaching material and CV pages.
- A modern research-oriented visual layer on top of the Academic Pages / Minimal Mistakes structure.

## Local Development

This project is compatible with GitHub Pages. For local preview on macOS, Ruby 3.1 is the most reliable option for this Jekyll stack.

```bash
brew install ruby@3.1
export PATH="/opt/homebrew/opt/ruby@3.1/bin:/opt/homebrew/lib/ruby/gems/3.1.0/bin:$PATH"
bundle config set path vendor/bundle
bundle install
bundle exec jekyll serve --config _config.yml,_config.dev.yml --host 127.0.0.1 --port 4000 --livereload
```

Then open:

```text
http://127.0.0.1:4000/
```

Useful localized URLs:

- English: `http://127.0.0.1:4000/`
- French: `http://127.0.0.1:4000/fr/`
- Japanese: `http://127.0.0.1:4000/ja/`

## Build Check

```bash
export PATH="/opt/homebrew/opt/ruby@3.1/bin:/opt/homebrew/lib/ruby/gems/3.1.0/bin:$PATH"
bundle exec jekyll build --config _config.yml,_config.dev.yml
```

## Troubleshooting

If the local page shows only unstyled HTML, the site is probably being served with production asset URLs or an old Jekyll process is still running.

```bash
lsof -ti tcp:4000 | xargs kill
export PATH="/opt/homebrew/opt/ruby@3.1/bin:/opt/homebrew/lib/ruby/gems/3.1.0/bin:$PATH"
bundle exec jekyll serve --config _config.yml,_config.dev.yml --host 127.0.0.1 --port 4000 --livereload
```

Then hard-refresh the browser.

## Content Structure

- `_pages/`: static pages and localized landing/navigation pages.
- `_publications/`: publication entries.
- `_talks/`: talk entries.
- `_teaching/`: teaching entries.
- `_data/navigation.yml`: navigation menus for English, French and Japanese.
- `_data/ui-text.yml`: translated UI labels.
- `_sass/_moonmodern.scss`: modern visual layer.
- `files/`: downloadable PDFs and other public files.

## Deployment

The site is deployed with GitHub Pages from the repository root. Changes merged into `main` trigger a GitHub Pages rebuild.

Local-only artifacts such as `_site/`, `.bundle/`, `.jekyll-cache/` and `vendor/bundle/` are ignored.
