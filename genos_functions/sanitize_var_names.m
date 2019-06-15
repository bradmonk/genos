function [VALID_VARS] = sanitize_var_names(VARS)
%{
sanitize_var_names(VARS)

  This function takes VARS, a list of strings or character vectors
  you want to use as TABLE.Properties.VariableNames. Built-in methods
  then proceed to sanitize the list to ensure variables are unique, 
  and valid MATLAB identifiers. The sanitization process flow...

1. Whitespace characters are removed. If a whitespace character is 
   followed by a lowercase letter, the letter is converted to the 
   corresponding uppercase character.

2. The character vectors are then checked to ensure only valid 
   MATLAB characters (A?Z, a?z, 0?9, _ )  are being used. Any
   invalid character found will be replaced by an underscore.

3. Determines if any duplicate strings exist, and appends an
   underscore and a number (e.g. 'var', 'var_1') to any duplicates.



```matlab

% EXAMPLE
%-------------

T = array2table(randi(6,6))

shitList = {'1var' '1var' ' 1var' 'my var' 'v#3' 'alpha and ?'}

validList = matlab.lang.makeValidName(shitList)

validListNoDupes = matlab.lang.makeUniqueStrings(validList)

T.Properties.VariableNames = validListNoDupes;

%-------------

```


%}



VALID_VARS = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(VARS));


end