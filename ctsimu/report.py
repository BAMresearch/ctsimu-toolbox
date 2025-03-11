# -*- coding: UTF-8 -*-
"""
PDF report generation.

Based on:
https://medium.com/@parveengoyal198/mastering-pdf-report-generation-with-reportlab-a-comprehensive-tutorial-part-2-c970ccd15fb6
"""

#from reportlab.pdfgen import canvas 
#from reportlab.pdfbase import pdfmetrics 
#from reportlab.lib import colors 
from datetime import datetime

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, NextPageTemplate, PageBreak, Table
from reportlab.pdfgen.canvas import Canvas
from reportlab.graphics.shapes import Drawing, Line, Rect, Ellipse
from reportlab.pdfbase import pdfmetrics
from reportlab.lib import colors 

class Report(BaseDocTemplate):
    """Template class for PDF report"""

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.bodyWidth = self.pagesize[0] - self.leftMargin - self.rightMargin
        self.bodyXCenter = self.leftMargin + .5 * self.bodyWidth
        self.headerPad = 6*mm
        self.header_x1 = self.leftMargin
        self.header_y1 = self.pagesize[1] - self.topMargin + self.headerPad
        self.footer_x1 = self.leftMargin
        self.footer_y2 = self.bottomMargin - self.headerPad

        self.header_font = 'Courier'
        self.header_fsize = 9
        h = pdfmetrics.TypeFace(self.header_font)
        self.header_fheight = 2 + ((h.ascent - h.descent) * self.header_fsize) / 1000
        self.headerLeft = ''
        self.headerCenter = ''
        self.headerRight = ''
        self.footerLeft = ''
        self.footerLeft2 = ''
        self.footerCenter = ''
        self.footerRight = ''
    
        # Define the page frames
        self.frame = Frame(
            self.leftMargin, self.bottomMargin, 
            self.bodyWidth, self.pagesize[1] - self.topMargin - self.bottomMargin,
            leftPadding=0, bottomPadding=0,
            rightPadding=0, topPadding=0,
            id='normal'
        )

        # Define the PageTemplate
        self.addPageTemplates([
            PageTemplate(
                id='normal',
                frames=[self.frame],
                onPage=self._header,
                onPageEnd=self._footer
            )
        ])

    def _header(self, canvas, doc):
        # Draw the header
        canvas.saveState() 
        canvas.setFont(self.header_font, self.header_fsize)
        yLine = self.header_y1
        canvas.line(self.header_x1, yLine, self.header_x1 + self.bodyWidth, yLine)
        yText = yLine + 2 * mm
        canvas.drawString(self.leftMargin, yText, self.headerLeft)
        canvas.drawCentredString(self.bodyXCenter, yText, self.headerCenter)
        canvas.drawRightString(self.leftMargin + self.bodyWidth, yText, self.headerRight)
        canvas.restoreState()

    def _footer(self, canvas, doc):
        # Draw the footer
        canvas.saveState() 
        canvas.setFont(self.header_font, self.header_fsize)
        yLine = self.footer_y2
        canvas.line(self.footer_x1, yLine, self.footer_x1 + self.bodyWidth, yLine)
        yText = yLine - 2 * mm - self.header_fheight
        canvas.drawString(self.footer_x1, yText, self.footerLeft)
        canvas.drawCentredString(self.footer_x1 + 0.5 * self.bodyWidth, yText, self.footerCenter)
        canvas.drawRightString(self.footer_x1 + self.bodyWidth, yText, self.footerRight)
        if canvas._pageNumber == 1:
            canvas.drawString(self.footer_x1, yText - self.header_fheight, self.footerLeft2)
        canvas.restoreState()
        # Use the right-string position for the page number in NumberedCanvas
        try:
            canvas.page_number_xpos = self.footer_x1 + self.bodyWidth
            canvas.page_number_ypos = yText
            canvas.page_number_font = self.header_font
            canvas.page_number_fsize = self.header_fsize
            canvas.page_number_text = 'Page %d of %d'
        except:
            pass

    def setHeader(self, left = '', center = '', right = ''):
        self.headerLeft = left
        self.headerCenter = center
        self.headerRight = right

    def setFooter(self, left = '', center = '', right = '', left2 = ''):
        self.footerLeft = left
        self.footerLeft2 = left2
        self.footerCenter = center
        self.footerRight = right

    def df2table(self, df):
        data = [df.columns.values.tolist()] + df.values.tolist()
        return Table(data, 45,
            repeatRows=1,
            style=[
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 7),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 4),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
                ('ALIGN', (0, 1), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
                ('FONTSIZE', (0, 1), (-1, -1), 9),
                ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)],
            # hAlign = 'LEFT',
            spaceBefore=10, spaceAfter=10)

class NumberedCanvas(Canvas):
    def __init__(self, *args, **kwargs):
        Canvas.__init__(self, *args, **kwargs)

        self._saved_page_states = []
        #print(self.getAvailableFonts())

        self.page_number_xpos = 190*mm
        self.page_number_ypos = 14*mm
        self.page_number_font = 'Helvetica'
        self.page_number_fsize = 9
        self.page_number_text = 'Page %d of %d'

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        """add page info to each page (page x of y)"""
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_number(num_pages)
            Canvas.showPage(self)
        Canvas.save(self)

    def draw_page_number(self, page_count):
        self.setFont(self.page_number_font, self.page_number_fsize)
        self.drawRightString(self.page_number_xpos, self.page_number_ypos,
            self.page_number_text % (self._pageNumber, page_count))


def main():
    # Create a new PDF document using the template
    pdf_doc = Report('example_page_template_header_footer.pdf', 
                    headerText = 'qwer',
                    pagesize=A4,
                    pageTemplates=[],
                    showBoundary=1,
                    leftMargin=27*mm,
                    rightMargin=20*mm,
                    topMargin=25*mm,
                    bottomMargin=25*mm,
                    allowSplitting=1,
                    title=None,
                    author=None,
                    _pageBreakQuick=1,
                    encrypt=None)

    pdf_doc.setHeader('l', 'c', 'r')
    pdf_doc.setFooter('l')

    # Add the content to the PDF document
    elements = [Paragraph('This is some content for the PDF document.'), PageBreak(), Paragraph('This is some content for the PDF document Page 2.')]

    pdf_doc.build(elements)

if __name__ == "__main__":
    main()