"""
Created on 10. feb. 2016

@author: skrinjar.luka@gmail.com

"""
import xlrd
import xlsxwriter

workbook = xlrd.open_workbook("test.xlsx")

#   get solution worksheet
worksheet = workbook.sheet_by_name("List1")



#    write data to excel
workbook = xlsxwriter.Workbook("test.xlsx")
#   add worksheet
worksheet = workbook.add_worksheet(self.__excel_worksheet)