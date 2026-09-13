# Research website design system

The Jekyll/Liquid architecture, public routes, scientific content and downloadable files remain the source of truth. This design replaces the accumulated moonmodern/refinements/research-design overrides with three component stylesheets.

## Tokens and typography

`_sass/_studio.scss` defines the shared palette, fonts, spacing, radii and motion. Use its custom properties in new components instead of introducing another global override stylesheet.

| Role | Token | Value |
| --- | --- | --- |
| Primary text | `--ink` | `#10233d` |
| Secondary text | `--muted` | `#52677f` |
| Action | `--blue` | `#155bd7` |
| Supporting accent | `--cyan` | `#07748a` |
| Canvas | `--canvas` | `#f4f7fb` |
| Surface | `--surface` | `#ffffff` |
| Divider | `--line` | `#dce5ef` |
| Panel corner | `--radius` | `16px` |

Manrope variable is self-hosted in Latin, Latin-extended and Greek subsets (49,400 bytes total), with system Japanese fallbacks and font-display: swap. The OFL license and source details are in assets/fonts/manrope. Base text is 16px, long desktop reading text 17px, article width 70ch. Major headings use a fluid scale, with larger typography reserved for the homepage. Decoration never conveys scientific information.

## Components

- `_sass/_studio.scss`: navigation, page grid, profile rail, home, footer, buttons, form states, responsive system and reduced-motion behavior.
- `_sass/_editorial.scss`: archive rows, scientific article headers, downloads, citations and previous/next navigation.
- `_sass/_profile-pages.scss`: shared CV and Contact pages.
- `_includes/masthead.html`: language-aware primary navigation; direct language links remain visible on small screens.
- `_includes/editorial-*.html`: shared article labels, resources and pagination.
- `_includes/profile-cv.html`: shared structure for the three localized CVs.

The main document is a grid above a footer in normal document flow. At 959px and below, the profile rail becomes a compact row and navigation becomes a disclosure panel. At 640px and below, the home art panel, research cards and CV sections rearrange for phone reading. JS-disabled navigation remains accessible. No section waits for an observer before becoming visible.

## Interactions and builds

`assets/js/site-enhancements.js` owns keyboard-accessible disclosures: native buttons, aria-expanded/aria-controls, Escape, outside-click dismissal and focus recovery. Motion respects prefers-reduced-motion. The visitor counter is unchanged.

`npm run build:js` uses Node built-ins to concatenate the preserved jQuery, FitVids and Magnific Popup sources plus the small initialization file. The legacy main.min.js route is retained for compatibility. No build-time package install is required. The obsolete greedy menu, Stickyfill and background footer timers have been removed from the served bundle. Original vendor files remain available in the repository. `node --check` verifies generated and standalone JavaScript.

Use the existing Ruby 3.1 environment and `bundle exec jekyll build --safe` for production. Do not edit generated HTML or PDF files to alter the design. All content downloads are preserved byte for byte.

## Validation and remaining content work

Review homepage, all three collection archives, article pages, all CV/contact languages, teaching resources and the 404 page at phone, tablet and desktop widths. Check navigation with keyboard, PDF/ZIP links, images, MathJax and the shared visitor counter. All source public routes are retained; the explicitly unpublished /cv-legacy/ remains unpublished.

The inherited /terms/ text describes analytics/comment services that are not enabled. Its wording warrants a separate content review; this design change does not introduce those services or alter that text.

## Futuristic atmosphere refinement

The full-width hero filaments are restored with a dark text veil, a cyan rim and a subtle portrait halo. The shared background uses static CSS halos and a grid masked toward the margins; the central reading area stays calm. Phone layouts simplify the background, and print removes all atmospheric pseudo-elements. No additional scripts, fonts or image requests are introduced. The numerical-linear-algebra illustration is now a decorative sparse matrix with column vectors, representing an abstract Ax = b without research data.

## Academic identity and illustrated portrait

The homepage name includes a responsive PhD suffix and the accessible name Ronan Dupont, PhD. The original avatar_2.png drawing is restored in the shared footer and in the localized Contact identity component. Its native size is 102 by 112 pixels; display sizes stay at or below that resolution. The footer image loads lazily and remains decorative beside the visible name. Contact uses a localized alternative description and keeps the existing contact details and professional links.
