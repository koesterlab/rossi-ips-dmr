import fitz

input_pdfs = snakemake.input.plots
output_pdf = snakemake.output[0]

# Größen der einzelnen PDFs
pages = []
for pdf in input_pdfs:
    doc = fitz.open(pdf)
    rect = doc[0].rect
    pages.append((doc, rect))

# Maximale Breite und Gesamthöhe
max_width = max(rect.width for _, rect in pages)
total_height = sum(rect.height for _, rect in pages)

out = fitz.open()
page = out.new_page(
    width=max_width,
    height=total_height,
)

y = 0

for doc, rect in pages:
    w = rect.width
    h = rect.height

    # Originalgröße beibehalten, linksbündig bei x=0
    page.show_pdf_page(
        fitz.Rect(0, y, w, y + h),
        doc,
        0,
    )

    y += h
    doc.close()

out.save(output_pdf)
out.close()
