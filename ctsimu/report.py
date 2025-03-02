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
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, NextPageTemplate, PageBreak
from reportlab.pdfgen.canvas import Canvas
from reportlab.graphics.shapes import Drawing, Line, Rect, Ellipse

class Report(BaseDocTemplate):
    """Template class for PDF report"""

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.headerLeft = ''
        self.headerCenter = ''
        self.headerRight = ''
        self.footerLeft = ''
        self.footerCenter = ''
        self.footerRight = ''
        self.boddyWidth = self.pagesize[0] - self.leftMargin - self.rightMargin
        self.boddyXCenter = self.leftMargin + .5 * self.boddyWidth
        self.headerPad = 6*mm

        # Define the page frames
        self.frame = Frame(
            self.leftMargin, self.bottomMargin, 
            self.boddyWidth, self.pagesize[1] - self.topMargin - self.bottomMargin,
            id='normal'
        )

        # Define the styles for header and footer
        self.styles = getSampleStyleSheet()
        self.header_style = self.styles['Code']
        self.footer_style = self.styles['Code']

        # Define the header and footer frames
        self.header_frame = Frame(
            self.leftMargin, self.pagesize[1] - self.topMargin + self.headerPad, self.boddyWidth, self.topMargin - self.headerPad - 6*mm,
            id='normal'
        )
        self.footer_frame = Frame(
            self.leftMargin, self.bottomMargin - self.headerPad - 6*mm, self.boddyWidth, 6*mm,
            id='footer'
        )

        # Define the PageTemplate
        self.addPageTemplates([
            PageTemplate(
                id='FirstPage',
                frames=[self.frame, self.header_frame, self.footer_frame],
                onPage=self._header,
                onPageEnd=self._footer
            )
        ])

    def _header(self, canvas, doc):
        # Draw the header
        self.header_style.alignment = 1  # center align the header text
        # header_text = Paragraph(f'My Header Text {self.leftMargin}', self.header_style)
        # header_text.wrapOn(canvas, self.header_frame.width, self.header_frame.height - 6*mm)
        # header_text.drawOn(canvas, self.header_frame.x1, self.header_frame.y1)
        yLine = self.header_frame.y1 - 1*mm
        canvas.line(self.header_frame.x1, yLine, self.header_frame.x1 + self.header_frame.width, yLine)
        yText = yLine + 6 * mm
        canvas.drawString(self.leftMargin, yText, self.headerLeft)
        canvas.drawCentredString(self.boddyXCenter, yText, self.headerCenter)
        canvas.drawRightString(self.leftMargin + self.boddyWidth, yText, self.headerRight)

    def _footer(self, canvas, doc):
        # Draw the footer
        self.footer_style.alignment = 2 #1  # center align the footer text
        footer_text = Paragraph("Page <seq id='PageNumber'/> of <seq id='TotalPages'/>", self.footer_style)
        footer_text.wrapOn(canvas, self.footer_frame.width, self.footer_frame.height)
        footer_text.drawOn(canvas, self.footer_frame.x1, self.footer_frame.y1)
        yLine = self.footer_frame.y1 + self.footer_frame.height + 1*mm
        canvas.line(self.footer_frame.x1, yLine, self.footer_frame.x1 + self.footer_frame.width, yLine)
        yText = yLine - 6 * mm
        canvas.drawString(self.footer_frame.x1, yText, self.footerLeft)
        canvas.drawCentredString(self.footer_frame.x1 + 0.5 * self.footer_frame.width, yText, self.footerCenter)
        canvas.drawRightString(self.footer_frame.x1 + self.footer_frame.width, yText, self.footerRight)

    def setHeader(self, left = '', center = '', right = ''):
        self.headerLeft = left
        self.headerCenter = center
        self.headerRight = right

    def setFooter(self, left = '', center = '', right = ''):
        self.footerLeft = left
        self.footerCenter = center
        self.footerRight = right

# # Create a new PDF document using the template
# pdf_doc = Report('example_page_template_header_footer.pdf', 
#                 headerText = 'qwer',
#                 pagesize=A4,
#                 pageTemplates=[],
#                 showBoundary=1,
#                 leftMargin=27*mm,
#                 rightMargin=20*mm,
#                 topMargin=25*mm,
#                 bottomMargin=25*mm,
#                 allowSplitting=1,
#                 title=None,
#                 author=None,
#                 _pageBreakQuick=1,
#                 encrypt=None)

# pdf_doc.setHeader('l', 'c', 'r')
# pdf_doc.setFooter('l')

# # Add the content to the PDF document
# elements = [Paragraph('This is some content for the PDF document.'), PageBreak(), Paragraph('This is some content for the PDF document Page 2.')]

# pdf_doc.build(elements)
