function [out1, out2] = exFun(inStr,inVar,inArr,inMat)
    % Short, one-line description. This is an example function.
    %{
    -----------------------------------------------------------------------
    DESCRIPTION:
        Longer description here in paragraph form. All text should fit
        within the standard MATLAB width unless the line includes a string.
        This is an example description. Lorem ipsum dolor sit amet, 
        consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
        labore et dolore magna aliqua. Ut enim ad minim veniam, quis 
        nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
        consequat. Duis aute irure dolor in reprehenderit in voluptate 
        velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint
        occaecat cupidatat non proident, sunt in culpa qui officia deserunt
        mollit anim id est laborum.
    -----------------------------------------------------------------------
    USAGE: [out1, out2] = funName(inStr,inVar,inArr,inMat)
        inStr:      str, string or character type input description
        inVar:      double, description of 1x1 variable (unit)
        inArr:      nx1 double, description of an array. Arrays shoulbe be
                    nx1, not 1xn whenever possible because MATLAB accesses
                    memory column-wise fastest
                        [var1 var2 var3]'
                        (unit unit unit)
        inMat:      nxm cell, description of matrix
                        {var1 var2 var3;
                         var4 var5 var6}
                        (unit unit unit;
                         unit unit unit)
    -----------------------------------------------------------------------
    OUTPUT: 
        out1:       type, description
        out2:       type, description
    -----------------------------------------------------------------------
    REFERENCES:
        - Citations here, preferrably IEEE format
    -----------------------------------------------------------------------
    NOTES:
        - Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do 
          eiusmod tempor incididunt ut labore et dolore magna aliqua.
        - Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do 
          eiusmod tempor incididunt ut labore et dolore magna aliqua.
    -----------------------------------------------------------------------
    AUTHOR: First Lastname
    -----------------------------------------------------------------------
    COPYTRIGHT: 2019 SLAB Group
    -----------------------------------------------------------------------
    TIMESTAMPS:
        - DD-MMM-YYYY: short description of change (initials)
    -----------------------------------------------------------------------
    %}
    
    %% Section1
    out1 = localFun(inVar);
    %% Section2
    out2 = inStr;

end

function out = localFun(in)
    % Short description
    % out = localFun(in)
    %{
    -----------------------------------------------------------------------
    INPUT:
        in:         type, description
    -----------------------------------------------------------------------
    OUTPUT:
        out:        type, description
    -----------------------------------------------------------------------
    NOTES:
        - Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do 
          eiusmod tempor incididunt ut labore et dolore magna aliqua.
    -----------------------------------------------------------------------
    %}
    
    out = in;

end