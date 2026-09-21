#!/usr/bin/env python3
"""SNAP amazon-meta.txt -> amazon-meta.tsv: one line per product with a title:
Id <TAB> group <TAB> salesrank <TAB> title <TAB> categories (full paths separated by ';', each path '|'-separated with
[ids]).  Usage: parse_amazon_meta.py raw/amazon-meta.txt amazon-meta.tsv"""
import sys

def main():
    src, dst = sys.argv[1], sys.argv[2]; n = 0
    with open(src, encoding='latin-1') as f, open(dst, 'w') as out:
        cur = None
        def flush():
            nonlocal n
            if cur and 'title' in cur:
                out.write(f"{cur['id']}\t{cur.get('group','')}\t{cur.get('salesrank','')}\t{cur['title']}\t{';'.join(cur.get('cats', []))}\n"); n += 1
        for line in f:
            if line.startswith('Id:'): flush(); cur = {'id': int(line.split()[1]), 'cats': []}
            elif cur is None: continue
            elif line.startswith('  title:'): cur['title'] = line[8:].strip().replace('\t', ' ')
            elif line.startswith('  group:'): cur['group'] = line[8:].strip()
            elif line.startswith('  salesrank:'): cur['salesrank'] = line[12:].strip()
            elif line.startswith('   |'): cur['cats'].append(line.strip())
        flush()
    print('products with a title', n)

if __name__ == '__main__':
    main()
