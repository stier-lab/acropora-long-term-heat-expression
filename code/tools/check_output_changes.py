#!/usr/bin/env python3
"""Sort uncommitted changes into real edits and rerun noise.

Run from the repo root after code/_run_all.R:

    python3 code/tools/check_output_changes.py            # report only
    python3 code/tools/check_output_changes.py --restore  # undo noise-only changes

"Noise" means a file whose only difference from the last commit is:
  - a PDF creation/modification timestamp or document ID
  - trailing whitespace or line endings in a text file
  - output/session_info.txt (it changes with whichever packages were loaded)
Everything else is listed as a real change to review and commit.
--restore runs `git checkout -- <file>` on noise-only files; nothing else is touched.
"""
import re, subprocess, sys

def git(*args, text=False):
    return subprocess.run(["git", *args], capture_output=True, text=text, check=False).stdout

def pdf_norm(b):
    b = re.sub(rb"/(CreationDate|ModDate)\s*\(D:[^)]*\)", b"", b)
    return re.sub(rb"/ID\s*\[\s*<[^>]*>\s*<[^>]*>\s*\]", b"", b)

def text_norm(b):
    return b"\n".join(l.rstrip() for l in b.replace(b"\r\n", b"\n").split(b"\n")).strip()

modified = [l[3:] for l in git("status", "--porcelain", text=True).splitlines() if l.startswith(" M")]
untracked = [l[3:] for l in git("status", "--porcelain", text=True).splitlines() if l.startswith("??")]
noise, real = [], []
for f in modified:
    old, new = git("show", f"HEAD:{f}"), open(f, "rb").read()
    ext = f.rsplit(".", 1)[-1].lower()
    if f == "output/session_info.txt":
        noise.append(f)
    elif ext == "pdf":
        (noise if pdf_norm(old) == pdf_norm(new) else real).append(f)
    elif ext in {"csv", "txt", "md", "r", "rmd", "html", "json", "tsv"}:
        (noise if text_norm(old) == text_norm(new) else real).append(f)
    else:
        real.append(f)

print(f"Real changes to review and commit ({len(real)}):")
for f in real: print("  ", f)
print(f"\nNoise only ({len(noise)}):")
for f in noise: print("  ", f)
if untracked:
    print(f"\nNew files not yet tracked ({len(untracked)}):")
    for f in untracked: print("  ", f)
if "--restore" in sys.argv and noise:
    subprocess.run(["git", "checkout", "--", *noise], check=True)
    print(f"\nRestored {len(noise)} noise-only files to the last commit.")
