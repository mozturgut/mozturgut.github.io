# AGENTS.md

## Cursor Cloud specific instructions

This is a static HTML/CSS/JS portfolio website (iPortfolio Bootstrap template). There is no build step, no package manager, no linting, and no automated tests.

### Running the dev server

Serve from the repository root with any static file server:

```
python3 -m http.server 8080
```

Then open `http://localhost:8080/` in Chrome.

### Key notes

- All vendor libraries (Bootstrap 5.3.3, AOS, Typed.js, PureCounter, Swiper, GLightbox, Isotope, etc.) are committed in `assets/vendor/` — no dependency installation is needed.
- The PHP contact form (`forms/contact.php`) requires a paid `php-email-form` library not included in the repo; it is non-functional.
- The R/Rmd files in the repo root are unrelated bioinformatics scripts, not part of the website.
