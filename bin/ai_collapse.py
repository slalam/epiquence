#!/usr/bin/env python3
"""Wrap MultiQC's global AI summary in a native <details> element so it collapses with no
JavaScript."""
import io
import re
import sys

REPORT = sys.argv[1] if len(sys.argv) > 1 else "multiqc_report.html"

DISCLAIMER = (
    '<div class="ai-disclaimer">For research and surveillance use only, not for '
    'clinical diagnosis. This summary is AI-generated and may be incomplete or '
    'inaccurate.</div>'
)

html = io.open(REPORT, encoding="utf-8", errors="ignore").read()
ORIGINAL = html


def searchable(s):
    """Return a copy with <style> and <script> bodies blanked out (same length), so that CSS/JS
    text can never be mistaken for real document markup."""
    def blank(m):
        return m.group(1) + (" " * len(m.group(2))) + m.group(3)
    s = re.sub(r"(<style\b[^>]*>)(.*?)(</style>)", blank, s, flags=re.S | re.I)
    s = re.sub(r"(<script\b[^>]*>)(.*?)(</script>)", blank, s, flags=re.S | re.I)
    return s


def strip_latex(s):
    r"""The model sometimes emits LaTeX ($\ge$, $\le$, $\beta$, \%) which MultiQC's markdown
    renderer does not process, so it shows up literally."""
    subs = [
        (r"\$\\ge\$", ">="), (r"\$\\geq\$", ">="),
        (r"\$\\le\$", "<="), (r"\$\\leq\$", "<="),
        (r"\$\\approx\$", "~"), (r"\$\\times\$", "x"),
        (r"\$\\alpha\$", "alpha"), (r"\$\\beta\$", "beta"),
        (r"\$\\gamma\$", "gamma"), (r"\$\\mu\$", "u"),
        (r"\\ge\b", ">="), (r"\\le\b", "<="),
        (r"\\%", "%"), (r"\\_", "_"),
    ]
    for pat, rep in subs:
        s = re.sub(pat, rep, s)
    # any leftover simple $...$ inline math: drop the delimiters
    s = re.sub(r"\$([^$<>]{1,40}?)\$", r"\1", s)
    return s


# de-LaTeX the AI summary response block(s) before wrapping
for _div in ("global_ai_summary_response", "global_ai_summary_detailed_analysis_response"):
    m_body = re.search(
        r'(<div[^>]*id="' + _div + r'"[^>]*>)(.*?)(</div>)', html, re.S
    )
    if m_body:
        cleaned = strip_latex(m_body.group(2))
        if cleaned != m_body.group(2):
            html = html[: m_body.start(2)] + cleaned + html[m_body.end(2):]

# The genuine wrapper is a <div ...id="global_ai_summary_wrapper"...> whose first child is <h5
# class="ai-summary-header">.
pat = re.compile(
    r'<div(\s[^>]*id="global_ai_summary_wrapper"[^>]*)>\s*'
    r'<h5[^>]*class="[^"]*ai-summary-header[^"]*"[^>]*>.*?</h5>',
    re.S,
)


def save_and_exit(msg):
    """Persist any changes made so far (e.g. the LaTeX cleanup) before exiting."""
    if html != ORIGINAL:
        io.open(REPORT, "w", encoding="utf-8").write(html)
        print(f"{msg} (LaTeX cleanup saved)")
    else:
        print(msg)
    sys.exit(0)


m = pat.search(searchable(html))
if not m:
    if '<summary class="ai-summary-header">' in searchable(html):
        save_and_exit("AI summary already wrapped")
    else:
        save_and_exit("AI summary wrapper not found")

attrs = m.group(1)
replacement = (
    f"<details open{attrs}>"
    '<summary class="ai-summary-header">'
    # Explicit caret, matching the outbreak/HGT accordions.
    '<span class="ai-caret">&#9654;</span>'
    'AI Interpretation of Results</summary>'
    + DISCLAIMER
)

# Find the </div> that closes this wrapper, by depth-counting from the match.
start = m.start()
body_start = m.end()
depth = 1          # we are inside the wrapper div
close_at = None
for tm in re.finditer(r"</?div\b[^>]*>", html[body_start:], re.I):
    if tm.group(0).lower().startswith("</div"):
        depth -= 1
        if depth == 0:
            close_at = body_start + tm.start()
            close_end = body_start + tm.end()
            break
    else:
        depth += 1

if close_at is None:
    save_and_exit("WARN: could not find closing tag for AI wrapper")

html = (
    html[:start]
    + replacement
    + html[body_start:close_at]
    + "</details>"
    + html[close_end:]
)

io.open(REPORT, "w", encoding="utf-8").write(html)
print("AI summary wrapped in <details> accordion")