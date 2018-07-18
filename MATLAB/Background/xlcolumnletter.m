function colLetter = xlcolumnletter(colNumber)
% Excel formats columns using letters.
% This function returns the letter combination that corresponds to a given
% column number.
% Limited to 702 columns
if( colNumber > 26*27 )
    error('XLCOLUMNLETTER: Requested column number is larger than 702. Need to revise method to work with 3 character columns');
else
    % Start with A-Z letters
    atoz        = char(65:90)';
      % Single character columns are first
      singleChar  = cellstr(atoz);
      % Calculate double character columns
      n           = (1:26)';
      indx        = allcomb(n,n);
      doubleChar  = cellstr(atoz(indx));
      % Concatenate
      xlLetters   = [singleChar;doubleChar];
      % Return requested column
      colLetter   = xlLetters{colNumber};
  end