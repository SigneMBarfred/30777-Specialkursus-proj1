# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:53:07 2023

@author: Oliver
"""

# Open the DAT file for reading
RAGERS = 'C:/Users/olive/Dropbox/5. semester/Special kursus/deboosted_completeness_sources.dat'
modified_RAGERS = 'RAGERS.dat'

# Attempt to read the file
with open(RAGERS, 'r') as input_file:
    lines = input_file.readlines()

# Initialize an empty list to store modified lines
modified_lines = []

# Process each line and multiply the values in the second column by 15
for line in lines: 
    columns = line.strip().split()
    if len(columns) >= 2:
        if columns[1].replace(".", "").isdigit():
            second_column_value = float(columns[1])
            modified_second_column = second_column_value * 15
            columns[1] = str(modified_second_column)
            modified_lines.append(' '.join(columns) + '\n')
        else:
            modified_lines.append(line)
    else:
        modified_lines.append(line)

# Write the modified lines to the output file
with open(modified_RAGERS, 'w') as output_file:
    output_file.writelines(modified_lines)




#Whats going on in each line
#for line in lines:: This starts a loop that iterates through each line in the lines list, representing each line of data in the input file.

#columns = line.strip().split(): This line splits the current line into a list of strings using whitespace as the separator. It first removes leading and trailing whitespace with strip() and then splits the line into words.

#if len(columns) >= 2:: This conditional statement checks if there are at least two columns in the columns list. If not, it indicates that there might be an issue with the line structure, and the line is kept unchanged in the modified_lines list.

#if columns[1].replace(".", "").isdigit():: This line checks if the second column (index 1 in the columns list) can be converted to a numeric value. It does this by removing any dots ('.') from the second column using replace(".", "") and then checking if the result consists only of digits using isdigit().

#second_column_value = float(columns[1]): If the second column is numeric, this line converts it to a float and assigns it to the second_column_value variable.

#modified_second_column = second_column_value * 15: This line multiplies the second_column_value by 15 and stores the result in the modified_second_column variable.

#columns[1] = str(modified_second_column): This line updates the second column in the columns list with the modified value, converting it back to a string.

#modified_lines.append(' '.join(columns) + '\n'): This line appends the modified line, created by joining the columns with spaces and adding a newline character, to the modified_lines list.

