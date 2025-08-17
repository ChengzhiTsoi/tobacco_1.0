# -*- coding: utf-8 -*-

def MFC(input_file:str, output_file:str, input_format="mol", output_format="cif"):
    from openbabel import pybel
    molecules = pybel.readfile(input_format, input_file)
    output_file_writer = pybel.Outputfile(output_format, output_file)
    for i, molecule in enumerate(molecules):
        output_file_writer.write(molecule)
    output_file_writer.close()