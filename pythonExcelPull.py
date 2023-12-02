
def pythonExcelPull(df, variable):
    found = False

    for column in df.columns:
        # Check if the variable name is in this column
        if df[column].eq(variable).any():
            # Get the index of the row where the variable name is found
            row_index = df.index[df[column] == variable].tolist()[0]

            # Assuming the value is in the next column
            next_column_index = df.columns.get_loc(column) + 1

            # Check if next column exists
            if next_column_index < len(df.columns):
                value_column = df.columns[next_column_index]
                return df.at[row_index, value_column]
               
    # If the variable name is not found
    return None
