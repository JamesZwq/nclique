#!/usr/bin/env python3
"""Detect long prose sentences in LaTeX paper sections.

Borrowed from skills/paper-architect/sentence-length-detector.md.
"""
import re, glob, os, sys

def strip_nested(text, cmd):
    pattern = r'\\' + cmd + r'\{'
    while True:
        m = re.search(pattern, text)
        if not m: break
        i = m.end()
        depth = 1
        while i < len(text) and depth > 0:
            if text[i] == '{': depth += 1
            elif text[i] == '}': depth -= 1
            i += 1
        text = text[:m.start()] + ' ' + text[i:]
    return text

def strip_for_analysis(text):
    text = re.sub(r'(?m)^\s*%.*$', '', text)
    text = re.sub(r'(?<!\\)%.*$', '', text, flags=re.MULTILINE)
    for env in ['algorithm', 'algorithm\\*', 'algorithmic',
                'table', 'table\\*', 'figure', 'figure\\*',
                'tabular', 'equation', 'equation\\*',
                'align', 'align\\*', 'gather', 'gather\\*',
                'enumerate', 'itemize', 'proof',
                'example', 'lemma', 'theorem', 'definition', 'corollary',
                'property', 'remark']:
        text = re.sub(r'\\begin\{' + env + r'\}.*?\\end\{' + env + r'\}',
                      ' BLOCK ', text, flags=re.DOTALL)
    text = re.sub(r'\\\[.*?\\\]', ' MATH ', text, flags=re.DOTALL)
    text = re.sub(r'\$\$.*?\$\$', ' MATH ', text, flags=re.DOTALL)
    text = re.sub(r'\$[^$]*\$', ' M ', text)
    for cmd in ['footnote', 'caption', 'Description']:
        text = strip_nested(text, cmd)
    text = re.sub(r'\\(cite|ref|label|cref|texorpdfstring|index|url|nocite)\{[^}]*\}',
                  ' ', text)
    text = re.sub(r'\\(section|subsection|subsubsection|stitle|sstitle|ssstitle|paragraph|chapter|expsection)\*?\{[^}]*\}',
                  ' SECTITLE ', text)
    text = re.sub(r'\\(emph|textbf|textit|texttt|underline|mathit|mathrm|mathbf|mathcal|mathsf|mbox)\{([^{}]*)\}',
                  r'\2', text)
    text = re.sub(r'\\[a-zA-Z]+\*?\{\s*\}', ' M ', text)
    text = re.sub(r'\\[a-zA-Z]+\*?', ' ', text)
    text = re.sub(r'[{}\\~^_]', ' ', text)
    return text

def find_long_sentences(directory, threshold):
    results = []
    for f in sorted(glob.glob(f'{directory}/*.tex')):
        # skip conflicted dropbox copies
        if 'conflicted' in f: continue
        with open(f) as fh:
            raw = fh.read()
        text = strip_for_analysis(raw)
        sentences = re.split(r'(?<=[.!?])\s+(?=[A-Z`a-z])', text)
        for s in sentences:
            s = re.sub(r'\s+', ' ', s).strip()
            if not s or 'BLOCK' in s or 'SECTITLE' in s:
                continue
            nw = len(s.split())
            if nw > threshold:
                results.append((f, nw, s[:300]))
    results.sort(key=lambda x: -x[1])
    return results

if __name__ == '__main__':
    threshold = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    directory = sys.argv[2] if len(sys.argv) > 2 else \
        '/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Sigmod2027Nuclear/sections'
    hits = find_long_sentences(directory, threshold)
    print(f'=== {len(hits)} prose sentences > {threshold} words ===\n')
    for f, nw, s in hits:
        print(f'[{nw}w] {os.path.basename(f)}: {s}')
        print()
