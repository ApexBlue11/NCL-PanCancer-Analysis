"""Build a pandoc reference.docx that renders tables legibly.

Pandoc's default reference document defines a Table style with no cell borders
(only a rule under the header row) and inherits body font size, so a wide
result table renders as unruled, oversized text that runs off the page. This
writes a reference document with bordered cells, a shaded header row, an 8 pt
table font and wider page margins, which is what the manuscript tables need.
"""
import os, re, shutil, subprocess, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "manuscript", "assets")
REF = os.path.join(ASSETS, "reference.docx")

TABLE_STYLE = """<w:style w:type="table" w:default="1" w:styleId="Table">
  <w:name w:val="Table"/>
  <w:basedOn w:val="TableNormal"/>
  <w:qFormat/>
  <w:pPr><w:spacing w:before="20" w:after="20" w:line="200" w:lineRule="auto"/></w:pPr>
  <w:rPr><w:sz w:val="16"/><w:szCs w:val="16"/></w:rPr>
  <w:tblPr>
    <w:tblInd w:w="0" w:type="dxa"/>
    <w:tblBorders>
      <w:top w:val="single" w:sz="6" w:space="0" w:color="666666"/>
      <w:left w:val="single" w:sz="6" w:space="0" w:color="AAAAAA"/>
      <w:bottom w:val="single" w:sz="6" w:space="0" w:color="666666"/>
      <w:right w:val="single" w:sz="6" w:space="0" w:color="AAAAAA"/>
      <w:insideH w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>
      <w:insideV w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>
    </w:tblBorders>
    <w:tblCellMar>
      <w:top w:w="40" w:type="dxa"/><w:left w:w="80" w:type="dxa"/>
      <w:bottom w:w="40" w:type="dxa"/><w:right w:w="80" w:type="dxa"/>
    </w:tblCellMar>
  </w:tblPr>
  <w:tblStylePr w:type="firstRow">
    <w:rPr><w:b/><w:sz w:val="16"/></w:rPr>
    <w:tcPr>
      <w:shd w:val="clear" w:color="auto" w:fill="EDEDED"/>
      <w:tcBorders><w:bottom w:val="single" w:sz="12" w:space="0" w:color="333333"/></w:tcBorders>
    </w:tcPr>
  </w:tblStylePr>
</w:style>"""


def main():
    os.makedirs(ASSETS, exist_ok=True)
    subprocess.run(["pandoc", "-o", REF, "--print-default-data-file", "reference.docx"],
                   check=True)

    tmp = REF + ".tmp"
    zin = zipfile.ZipFile(REF)
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/styles.xml":
                s = data.decode("utf-8")
                s = re.sub(r'<w:style w:type="table" w:default="1" w:styleId="Table">.*?</w:style>',
                           TABLE_STYLE, s, count=1, flags=re.S)
                data = s.encode("utf-8")
            elif item.filename == "word/document.xml":
                # Widen the page so wide tables fit. The default reference
                # document has a bare <w:sectPr> with no page size or margins,
                # so these must be inserted rather than substituted.
                s = data.decode("utf-8")
                page = ('<w:pgSz w:w="12240" w:h="15840"/>'
                        '<w:pgMar w:top="907" w:right="907" w:bottom="907" '
                        'w:left="907" w:header="708" w:footer="708" w:gutter="0"/>')
                if "<w:pgMar" in s:
                    s = re.sub(r'<w:pgSz[^/]*/>\s*', '', s)
                    s = re.sub(r'<w:pgMar[^/]*/>', page, s)
                else:
                    s = s.replace("<w:sectPr>", "<w:sectPr>" + page, 1)
                data = s.encode("utf-8")
            zout.writestr(item, data)
    zin.close()
    shutil.move(tmp, REF)

    z = zipfile.ZipFile(REF)
    s = z.read("word/styles.xml").decode("utf-8")
    assert "tblBorders" in s, "table borders were not injected"
    print("wrote %s" % REF)
    print("  tblBorders present : %s" % ("tblBorders" in s))
    print("  header shading     : %s" % ("EDEDED" in s))
    print("  table font 8pt     : %s" % ('<w:sz w:val="16"/>' in s))


if __name__ == "__main__":
    main()
