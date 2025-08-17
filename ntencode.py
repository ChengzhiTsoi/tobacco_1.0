# -*- coding: utf-8 -*-

def deldot(filename, dotname):
    '''

    :param filename: DESCRIPTION
    :type filename: TYPE
    :param dotname: DESCRIPTION
    :type dotname: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    '''
    import os
    name=[]
    for i in filename:
        portion = os.path.splitext(i)
        if portion[1] == dotname:
            name.append(portion[0])
    name.sort()
    return name

def encode(sub_name,filename):
    import operator
    code=[]
    for i in filename:
        code_raw = []
        for j in sub_name:
            if operator.contains(i,j)==True:
                m=1
                code_raw.append(m)           
            else:
                m=0
                code_raw.append(m) 
        code.append(code_raw)
    return code