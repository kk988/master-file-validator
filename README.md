# master-file-validator
Script that will make sure the master file is correctly formatted.

## How to run this:

`Rscript Berger_Format_Validator.R "design.file='test_files/test_file_Multi_errors.txt'"`
`Rscript Berger_Format_Validator.R "design.file='test_files/test_file_CORRECT.txt'"`

You can do full path or not full path.

If the script finds any basic errors, such as not having enough columns, or
not having numbers in a column that needs to have numbers for downstream
validations, the script will exit.

If there are no basic errors, the script will then check for errors in general
formatting and logic checking. They will be printed out together, since one of
these errors will not affect any downstream validations.

