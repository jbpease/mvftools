python3 makeargparse.py
cat intro.md prog-desc.md footer.md > mvftools.md
pandoc -f markdown -t html -o mvftools.html --standalone --toc --toc-depth=2 mvftools.md
pandoc -f markdown -t latex -o mvftools.pdf --standalone --toc --toc-depth=2 mvftools.md
