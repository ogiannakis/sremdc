;+
; NAME:
;        STRREPLACE
;
; PURPOSE:
;        The STRREPLACE procedure replaces the contents of one string
;        with another.  The first occurrence of the search substring, Find
;        within the source string, String is replaced by the string,
;        Replacement.
;
; CATEGORY:
;        String Processing.
;
; CALLING SEQUENCE:
;
;        STRREPLACE, String, Find, Replacement
;
; INPUTS:
;        String:   The string to have substring(s) replaced.  If String is
;                  an array, Find is replaced by Replacement in the first
;                  occurrence of Find of every element of the array.
;
;        Find:     The scalar substring to be replaced. If this argument is
;                  not a string, it is converted using IDL's default
;                  formatting rules.
;
;        Replacement:   A scalar string to replace the Find substring. If
;                  this argument is not a string, it is converted using IDL's
;                  default formattting rules.
;
; EXAMPLE:
;
;        If the variable A contains the string "IBM is fun", the
;        substring "IBM" can be replaced with the string "Microsoft"
;        by entering:
;
;        STRREPLACE, A, 'IBM', 'Microsoft'
;
; MODIFICATION HISTORY:
;        Written by:    Han Wen, June 1995.
;-
function STRREPLACE, Strings, Find1, Replacement1

;   Check integrity of input parameter

Strings2=Strings

         NP        = N_PARAMS()
         if (NP ne 3) then message,'Must be called with 3 parameters, '+$
                   'Strings, Find, Replacement'

         sz        = SIZE(Strings)
         ns        = n_elements(sz)
         if (sz(ns-2) ne 7) then message,'Parameter must be of string type.'

         Find      = STRING(Find1)
         pos       = STRPOS(Strings,Find)
         here      = WHERE(pos ne -1, nreplace)

         if (nreplace eq 0) then stop

         Replacement=STRING(Replacement1)
         Flen      = strlen(Find)
         for i=0,nreplace-1 do begin

              j         = here(i)
              prefix    = STRMID(Strings(j),0,pos(j))
              suffix    = STRMID(Strings(j),pos(j)+Flen,$
                                       strlen(Strings(j))-(pos(j)+Flen))
              Strings2(j) = prefix + replacement + suffix
         endfor

exit:

return,strings2

end
;---------------------------------------------------;
;---------------SUPPORTING FUNCTIONS----------------;
;-----------------------START-----------------------;
;---------------------------------------------------;
Function Type, x

;+
; NAME:
;	TYPE
; VERSION:
;	3.0	
; PURPOSE:
;	Finds the type class of a variable.
; CATEGORY:
;	Programming.
; CALLING SEQUENCE:
;	Result = TYPE(X)
; INPUTS:
;    X
;	Arbitrary, doesn't even need to be defined.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;	None.
; OUTPUTS:
;	Returns the type of X as a long integer, in the (0,9) range.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Extracts information from the SIZE function.
; MODIFICATION HISTORY:
;	Created 15-JUL-1991 by Mati Meron.
;-

    dum = size(x)
    return, dum(dum(0) + 1)
end
;------------------------------------------------------------------

Function Isnum, x, double = doub, complex = comp, type = typ

;+
; NAME:
;	ISNUM
; VERSION:
;	3.0
; PURPOSE:
;	Checks whether the input is a number.
; CATEGORY:
;	Programming.
; CALLING SEQUENCE:
;	Result = ISNUM(X)
; INPUTS:
;    X
;	Arbitrary, doesn't even have to exist.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    /DOUBLE
;	Switch.  If set the result is 1 only if X is DOUBLE or DCOMPLEX.
;    /COMPLEX
;	Switch.  If set the result is 1 only if X is COMPLEX or DCOMPLEX.
;    TYPE
;	Optional output.  See below.
; OUTPUTS:
;	Returns 1 if X is number, 0 otherwise.  Output type is byte.
; OPTIONAL OUTPUT PARAMETERS:
;    TYPE
;	The name of the variable to receive the numeric code of the type of X.
;	Included for convenience to save an additional call to TYPE.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Straightforward.  Using TYPE from MIDL.
; MODIFICATION HISTORY:
;	Created 15-JUN-1995 by Mati Meron.
;	Modified 5-MAY-1996 by Mati Meron.  Added keywords DOUBLE, COMPLEX and
;	TYPE.
;-

    numtyps = [1,2,3,4,5,6,9]
    typ = Type(x)
    res = (where(numtyps eq typ))(0) ge 0
    if keyword_set(doub) then res = res and (typ eq 5 or typ eq 9)
    if keyword_set(comp) then res = res and (typ eq 6 or typ eq 9)

    return, res
end
;-------------------------------------------------------------------------------
Function Arreq, arr1, arr2, warn = wn, novalue = nov

;+
; NAME:
;	ARREQ
; VERSION:
;	3.0
; PURPOSE:
;	Compares arrays for equality.  The arrays qualify as equal if:
;	    1) They are of the same general type (num., char., or struct.).
;	    2) Number of dimensions is the same.
;	    3) Size of each dimension is the same.
;	    4) Respective elements are equal.
; CATEGORY:
;	Mathematical Function (general).
; CALLING SEQUENCE:
;	Result = ARREQ( ARR1, ARR2 [, keywords])
; INPUTS:
;    ARR1, ARR2
;	Arrays, type and number of dimensions arbitrary.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    /WARN
;	Switch. If set, a warning message is issued for incompatible data types.
;    /NOVALUE
;	Switch.  If set, only number of elements and structure are compared.
; OUTPUTS:
;	Returns 1 if the arrays are equal, 0 otherwise.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Uses the SIZE function and ISNUM from MIDL to obtain information about 
;	the arrays.  Compares, in order, number of dimensions, size of each 
;	dimension, data types, and (unless NOVALUE is set) individual elements.
; MODIFICATION HISTORY:
;	Created 15-JUL-1991 by Mati Meron.
;	Modified 30-AUG-1998 by Mati Meron.  Scope broadened to pointer arrays.
;-

    fsiz = size(arr1)
    ssiz = size(arr2)
    if fsiz(0) eq ssiz(0) then ndim = fsiz(0) else return, 0b
    for i = 1, ndim do if fsiz(i) ne ssiz(i) then return, 0b
    fnum = Isnum(arr1, type = ftyp)
    snum = Isnum(arr2, type = styp)
    if not ((fnum and snum) or (ftyp eq styp)) then begin
	if keyword_set(wn) then message, 'Incompatible data types!', /continue
	return, 0b
    endif

    if keyword_set(nov) then return, 1b else return, min(arr1 eq arr2)
end
;-------------------------------------------------------------------------------
Function Cast, x, low, high, fix = fix

;+
; NAME:
;	CAST
; VERSION:
;	3.0
; PURPOSE:
;	Generalized type casting.  Converts variables whose type code is out 
;	of the range [LOW,HIGH] into this range.
; CATEGORY:
;	Programming (type conversion).
; CALLING SEQUENCE:
;	Result = CAST( X, [LOW [,HIGH]])
; INPUTS:
;    X
;	Numerical, arbitrary, or a character representation of a number(s).
;    LOW
;	Number representing a type code, range (1:9).  If greater than 9, it is
;	set to 9.  If less then 1, or not given, it is set to 1.
; OPTIONAL INPUT PARAMETERS:
;    HIGH
;	Type code, same as LOW.  Default value is 9.  If provided and less then
;	LOW, it is set to LOW.
; KEYWORD PARAMETERS:
;    /FIX
;	Switch.  If set, the output is filtered through FPU_FIX, eliminating
;	floating underflow errors.
; OUTPUTS:
;	If the type of X is < LOW, CAST returns X converted to type LOW.
;	If the type of X is > HIGH, CAST returns X converted to type HIGH.
;	Otherwise CAST returns X.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	1)  An attempt to convert a string which is NOT a character 
;	    representation of a number into a numeric type will yield error.
;	2)  X cannot be a structure or pointer, but can be a structure element.
;	3)  The value 8 for either LOW or HIGH is not allowed (since it 
;	    corresponds to structure type).  Value of 10 and above is ignored.
; PROCEDURE:
;	Identifies the type of X, and if out of the range given by [LOW,HIGH]
;	calls the proper conversion routine using the system routine 
;	CALL_FUNCTION.  Also uses FPU_FIX and ISNUM from MIDL.
; MODIFICATION HISTORY:
;	Created 25-DEC-1991 by Mati Meron.
;	Modified 15-JUN-1995 by Mati Meron to accept the new DOUBLECOMPLEX type.
;	Modified 25-SEP-1998 by Mati Meron.  Underflow filtering added.
;-

    on_error, 1
    conv = ['nada', 'byte', 'fix', 'long', 'float', 'double', 'complex', $
	    'string', 'nonap', 'dcomplex']
    if n_elements(low)  eq 0 then ilo = 1 else ilo = 1   > fix(low)  < 9
    if n_elements(high) eq 0 then ihi = 9 else ihi = ilo > fix(high) < 9

    inum = Isnum(x, type = ityp) or ityp eq 7
    if ilo eq 8 or ihi eq 8 or not inum then message, "Can't do that!" else $
    if ityp lt ilo then res = call_function(conv(ilo),x) else $
    if ityp gt ihi then res = call_function(conv(ihi),x) else res = x

    if keyword_set(fix) then return, FPU_fix(res) else return, res
end
;-----------------------------------------------------------------------------
Function Default, x, y, strict = strit, dtype = deft, low = lot, high = hit

;+
; NAME:
;	DEFAULT
; VERSION:
;	3.0
; PURPOSE:
;	Provides an automatic default value for nondefined parameters.
; CATEGORY:
;	Programming.
; CALLING SEQUENCE:
;	Result = DEFAULT( X, Y [, keywords])
; INPUTS:
;    X, Y
;	Arbitrary, at least one needs to be defined.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    /STRICT
;	Switch.  If set, X is considered defined only if it is of the same type 
;	as Y.
;    /DTYPE
;	Switch.  If set, the result will be typecast into the type of Y.  
;	Explicit settings for LOW and/or HIGH (see below) override DTYPE.
;    LOW
;	Numeric value between 1 to 9 (8 is excluded). If given, the result is 
;	of type >= LOW.
;    HIGH
;	Numeric value between 1 to 9 (8 is excluded). If given, the result is 
;	of type <= HIGH.
; OUTPUTS:
;	X if it is defined, otherwise Y.  
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	All type casting is bypassed if the result is of type 8 (STRUCTURE) or
;	higher than 9.
; PROCEDURE:
;	Uses the functions CAST, ISNUM and TYPE from MIDL.
; MODIFICATION HISTORY:
;	Created 15-JUL-1991 by Mati Meron.
;	Modified 15-NOV-1993 by Mati Meron.  The keyword TYPE has been replaced
;	by STRICT.  Added keywords DTYPE, LOW and HIGH.
;-

    on_error, 1
    xtyp = Type(x)
    ytyp = Type(y)

    if not (xtyp eq 0 or keyword_set(strit)) then atyp = xtyp else $
    if ytyp ne 0 then atyp = ytyp else message,'Insufficient data!'

    if xtyp eq atyp then res = x else res = y

    if keyword_set(deft) then begin
	if n_elements(lot) eq 0 then lot = ytyp
	if n_elements(hit) eq 0 then hit = ytyp
    end

    if Isnum(res) or atyp eq 7 then return, Cast(res,lot,hit) else return, res
end
;----------------------------------------------------------------------------------
Function One_of, v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, value = val

;+
; NAME:
;	ONE_OF
; VERSION:
;	3.0
; PURPOSE:
;	Called with up to 8 variables V_0 through V_7 , ONE_OF checks which 
;	variable is defined (only one is supposed to be).
; CATEGORY:
;	Programming.
; CALLING SEQUENCE:
;	Result = ONE_OF( V_0 [,V_1, ... V_7] [, VALUE = VAL])
; INPUTS:
;    V_0 through V_7
;	Arbitrary, all are optional.
; OPTIONAL INPUT PARAMETERS:
;	See above.
; KEYWORD PARAMETERS:
;    VALUE
;	Optional output, see below.
; OUTPUTS:
;	Returns the serial number (range 0 through 7) of the defined variable,
;	or -1 if none is defined.  If more than one variable is defined, ONE_OF
;	issues an error message and returns to the main level.
; OPTIONAL OUTPUT PARAMETERS:
;    VALUE
;	The name of the variable to receive the value of the single defined
;	variable, or a null string if none is defined.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	Currently ONE_OF is restricted to a maximum of 8 variables.  If needed,
;	the number can be increased.
; PROCEDURE:
;	Straightforward.
; MODIFICATION HISTORY:
;	Created 15-JUL-1991 by Mati Meron.
;	Modified 30-JUL-1991 by Mati Meron.  The dependence of the original 
;	code on the EXECUTE system routine has been eliminated in order to 
;	assure compatibility with the OUTPUT routine.
;	Modified 15-NOV-1993 by Mati Meron.  Since IDL now allows for recursive
;	calls to EXECUTE, the original code has been restored.
;-

    on_error, 1
    vnams = ['v_0','v_1','v_2','v_3','v_4','v_5','v_6','v_7']
    exlist = lonarr(8)
    exind = -1l
    val = ''
    
    for i = 0, n_params() - 1 do idum = $
	execute('exlist(i) = n_elements(' + vnams(i) + ')') 
    wex = where(exlist gt 0, nex)
    if nex eq 1 then begin
	 exind = wex(0)
	idum = execute('val = ' + vnams(exind))
    endif else if nex gt 1 then message, 'Only one variable may be defined!'
    return, exind
end
;----------------------------------------------------------------------------------
Function Extrema, x, min_only = mino, max_only = maxo, ceiling = ceil, $
    threshold = tre, signature = sig, number = num

;+
; NAME:
; EXTREMA
; VERSION:
; 3.0
; PURPOSE:
; Finding all local minima and maxima in a vector.
; CATEGORY:
; Mathematical Function (array).
; CALLING SEQUENCE:
; Result = EXTREMA( X [, keywords])
; INPUTS:
;    X
; Numerical vector, at least three elements.
; OPTIONAL INPUT PARAMETERS:
; None.
; KEYWORD PARAMETERS:
;    /MIN_ONLY
; Switch.  If set, EXTREMA finds only local minima.
;    /MAX_ONLY
; Switch.  If set, EXTREMA finds only local maxima.
;    THRESHOLD
; A nonnegative value.  If provided, entries which differ by less then
; THRESHOLD are considered equal.  Default value is 0.
;    /CEILING
; Switch.  Determines how results for extended extrema (few consecutive 
; elements with the same value) are returned.  See explanation in OUTPUTS.
;    SIGNATURE
; Optional output, see below.
;    NUMBER
; Optional output, see below.
; OUTPUTS:
; Returns the indices of the elements corresponding to local maxima 
; and/or minima.  If no extrema are found returns -1.  In case of 
; extended extrema returns midpoint index.  For example, if 
; X = [3,7,7,7,4,2,2,5] then EXTREMA(X) = [2,5].  Note that for the 
; second extremum the result was rounded downwards since (5 + 6)/2 = 5 in
; integer division.  This can be changed using the keyword CEILING which 
; forces upward rounding, i.e. EXTREMA(X, /CEILING) = [2,6] for X above.
; OPTIONAL OUTPUT PARAMETERS:
;    SIGNATURE
; The name of the variable to receive the signature of the extrema, i.e.
; +1 for each maximum and -1 for each minimum.
;    NUMBER
; The name of the variable to receive the number of extrema found.  Note
; that if MIN_ONLY or MAX_ONLY is set, only the minima or maxima, 
; respectively, are counted.
; COMMON BLOCKS:
; None.
; SIDE EFFECTS:
; None.
; RESTRICTIONS:
; None.
; PROCEDURE:
; Straightforward.  Calls ARREQ, DEFAULT and ONE_OF from MIDL.
; MODIFICATION HISTORY:
; Created 15-FEB-1995 by Mati Meron.
; Modified 15-APR-1995 by Mati Meron.  Added keyword THRESHOLD.
;-

    on_error, 1
    siz = size(x)
    if siz(0) ne 1 then message, 'X must be a vector!' else $
    if siz(1) lt 3 then message, 'At least 3 elements are needed!'

    len = siz(1)
    res = replicate(0l,len)
    sig = res
    both = One_of(mino,maxo) eq -1
    cef = keyword_set(ceil)
    tre = Default(tre,0.,/dtype) > 0

    xn = [0, x(1:*) - x(0:len-2)]
    if tre gt 0 then begin
  tem = where(abs(xn) lt tre, ntem)
  if ntem gt 0 then xn(tem) = 0
    endif
    xp = shift(xn,-1)
    xn = xn(1:len-2)
    xp = xp(1:len-2)

    if keyword_set(mino) or both then begin
  fir = where(xn lt 0 and xp ge 0, nfir) + 1
  sec = where(xn le 0 and xp gt 0, nsec) + 1
  if nfir gt 0 and Arreq(fir,sec) then begin
      res(fir) = fir
      sig(fir) = -1
  endif else begin
      if nfir le nsec then begin
    for i = 0l, nfir-1 do begin
        j = (where(sec ge fir(i)))(0)
        if j ne -1 then begin
      ind = (fir(i) + sec(j) + cef)/2
      res(ind) = ind
      sig(ind) = -1
        endif
    endfor
      endif else begin
    for i = 0l, nsec-1 do begin
        j = (where(fir le sec(i), nj))((nj-1) > 0)
        if j ne -1 then begin
      ind = (sec(i) + fir(j) + cef)/2
      res(ind) = ind
      sig(ind) = -1
        endif
    endfor
      endelse
  endelse
    endif

    if keyword_set(maxo) or both then begin
  fir = where(xn gt 0 and xp le 0, nfir) + 1
  sec = where(xn ge 0 and xp lt 0, nsec) + 1
  if nfir gt 0 and Arreq(fir,sec) then begin
      res(fir) = fir
      sig(fir) = 1
  endif else begin
      if nfir le nsec then begin
    for i = 0l, nfir-1 do begin
        j = (where(sec ge fir(i)))(0)
        if j ne -1 then begin
      ind = (fir(i) + sec(j) + cef)/2
      res(ind) = ind
      sig(ind) = 1
        endif
    endfor
      endif else begin
    for i = 0l, nsec-1 do begin
        j = (where(fir le sec(i), nj))((nj-1) > 0)
        if j ne -1 then begin
      ind = (sec(i) + fir(j) + cef)/2
      res(ind) = ind
      sig(ind) = 1
        endif
    endfor
      endelse
  endelse
    endif

    res = where(res gt 0, num)
    sig = sig(res > 0)

    return, res
end
;-------------------------------------------------------------------------------
function Fsort, Array, Asort, INFO=info, DESCENDING=descend
;+
; NAME:
; FSORT
; PURPOSE:
; Function to sort data into ascending order,
; original subscript order is maintained when values are equal (FIFO).
; (unlike IDL sort routine alone,
;     which may rearrange the order of equal values)
; CALLING SEQUENCE:  
; result = fsort( array, [asort] )
; INPUT:
; Array - array to be sorted
; KEYWORDS:
; /INFO = optional keyword to cause brief message about # equal values.
; /DESCENDING = descending order instead of the default ascending order.
; OUTPUT:
; Subscripts which give sorted order are returned as function value.
; OPTIONAL OUTPUT:
; Asort - sorted array
; METHOD:
; uses WHERE to find equal clumps, instead of looping with IF ( EQ ).
; HISTORY:
; written by F. Varosi NASA/GSFC 1990.
;-
  N = N_elements( Array )
  if (N EQ 1) then return,[0]    ;Only 1 element

  if (N LE 0) then begin
    message,"expecting an array as input",/INFO
    retall
     endif

  subs = sort( Array )
  Asort = Array(subs)

; now sort subscripts into ascending order
; whenever clumps of equality are found in Asort:

  weq = where( shift( Asort, -1 ) EQ Asort, Neq )

  if keyword_set( info ) then message, strtrim( Neq,2 ) + $
            " equal values Located",/INFO

  if (Neq GE N) then begin

    subs = Lindgen( N )
    Asort = Array

   endif else if (Neq GT 0) then begin

    if (Neq GT 1) then begin        ;find clumps of equality

      wclump = where( (shift( weq, -1 ) - weq) GT 1, Nclump )
      Nclump = Nclump + 1

      endif else Nclump = 1

    if (Nclump LE 1) then begin
      Clump_Beg = 0
      Clump_End = Neq-1
      endif else begin
      Clump_Beg = [0,wclump+1]
      Clump_End = [wclump,Neq-1]
       endelse

    weq_Beg = weq( Clump_Beg )        ;subscript ranges
    weq_End = weq( Clump_End ) + 1    ; of Asort equalities.

    if keyword_set( info ) then message, strtrim( Nclump, 2 ) + $
        " clumps of equal values Located",/INFO

    for ic = 0L, Nclump-1 do begin    ;sort each clump.

      subic = subs( weq_Beg(ic) : weq_End(ic) )
            subs( weq_Beg(ic) ) = subic( sort( subic ) )
      endfor

    if keyword_set( descend ) then subs = rotate( subs, 2 )
    if N_params() GE 2 then Asort = Array(subs) ;resort array.
     endif

return, subs
end
;-------------------------------------------------------------------------------
;+
; NAME:
; InterLeave
; PURPOSE:
; Determine the interleaving of two arrays of numbers,
; by returning subscripts of the interleaving elements,
; basically binning the array B into bins defined by array A.
; Note that the input arrays can be in any order (non-sorted), but
; for using result of this function you probably want A sorted first.
; CALLING:
; Locs = InterLeave( A, B )
; INPUTS:
; A = the reference grid array of numbers.
; B = array of numbers to interleave with A.
; OUTPUT:
; Function returns an integer array of same length as input B,
; containing for each element of B the subscript of the element in A
; which is less than (or equal) and closest to that element of B.
; If there is no element of A less than or equal to an element of B
; that subscript is set to -1.
; KEYWORD:
; NBINS = optional output, # how many bins of grid A were occupied.
; EXTERNAL CALLS:
; function Fsort  (to sort and keep original order if equal).
; PROCEDURE:
; Proceed as in the beginning of pro match (by D.Lindler),
; but then find interleave locations by observing jumps in flag array.
; HISTORY:
; Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
;-

function InterLeave, a, b, NBINS=nstart

  na = N_elements( a )
  nb = N_elements( b )

  if (na LE 0) OR (nb LE 0) then begin
    message,"must supply two arrays for interleaving",/INFO
    print,"syntax:  B_in_A = interleave( A, B )"
    return,[0]
     endif

  if min( B ) LT min( A ) then begin
    if (na GT 32768) then Locs = Lonarr( nb )  $
         else Locs = intarr( nb )
    Locs(*) = -1
    w = where( B GE min( A ), nw )
    if (nw GT 0) then Locs(w) = InterLeave( A, B(w), NBINS=nstart )
    return, Locs
     endif

  sub = Fsort( [a,b] )      ;combine a and b and sort.
  ind = [ Lindgen( na ), Lindgen( nb ) ]  ;combined list of indices.
  ind = ind(sub)        ;same sort order

  flag = [ bytarr( na ), replicate( 1B, nb ) ]  ;to indicate which array
  flag = flag(sub)        ;same sort order.
  sub = 0
  flags1 = shift( flag, -1 )

; Find interleave locations by observing jumps in flag array:

  wstart = where( flags1 GT flag, nstart )
  wstop = where( flags1 LT flag, nstop )

  if (nstart NE nstop) then begin
    message,"possible problem:   nstart NE nstop",/INFO
    print, nstart, nstop, string( 7B )
     endif

  flag = 0
  flags1 = 0

  if (na GT 32768) then  Locs = Lonarr( nb )  else  Locs = intarr( nb )
  ws1 = wstart + 1
  for i=0,nstop-1 do Locs(ind(ws1(i):wstop(i))) = ind(wstart(i))

return, Locs
end
;--------------------------------------------------------------------------------
;+
; NAME:
; FINTERPOL
; PURPOSE:
; Linear interpolation of a table of function values onto a new grid.
; Finterpol is faster than the old IDL lib routine interpol
; when the arrays (xold & xnew) have more than about 10 elements,
; and speedup factor increases to > 15 when number of elements > 100.
; This is acheived by first calling function InterLeave( xold, xnew ),
; which then allows subsequent computations to be vectorized.
; Intermediate computations are stored in a common block so that
; after the first call, further interpolations of different functions
; evaluated on the same grids can be performed even faster (10 X)
; by omitting the xold & xnew arguments in further calls. The array
; of points at which function is known must be strictly increasing,
; but the new interpolation points can be in any order. If the new
; grid points are beyond the range of old points extrapolation is
; performed, and a message is printed (unless /QUIET is set).
; See keywords for options to interpolate exponential or power-law
; functional behaviours.
;
; CALLING:
; fxnew = finterpol( fxold, xold, xnew )
;
; INPUTS:
; fxold = array of function values at the points xold.
;
; xold = array of points (strictly increasing) where func is evaluated.
;
; xnew = array of new points at which linear interpolations are desired.
;
; KEYWORDS:
; /QUIET : do not check or print any messages if extrapolation is done.
;   Default is to give warning messages about # points extrapolated.
;
; /INITIALIZE : use only the xold & xnew arrays to compute arrays
;   that are kept in the common block for use in fast repetitive
;   interpolations of different functions from/to same grids.
;   No interpolation is performed and just status (0/1) is returned.
;
; /EXPONENTIAL_INTERP: perform interpolation/extrapolation in Log-Linear
;   space, thereby giving correct result for exponential function.
;   (Function values are then assumed > 1e-37).
;
; /POWER_LAW_INTERP: perform interpolation/extrapolation in Log-Log space,
;   thereby giving correct result for power-law function.
;   (Function and X coordinate values are then assumed > 1e-37).
;
; OUTPUTS:
; Function returns array of linear interpolates at new grid points,
; or it returns just the number 1 if /INIT is set, or 0 if calling error.
; COMMON BLOCKS:
; common finterpol
; common finterpol2
; EXTERNAL CALLS:
; function InterLeave
; function scalar
; PROCEDURE:
; Call function InterLeave to find old grid points which bracket
; the new grid points, and then interpolate in vectorized form.
; Take natural log of function values for exponential interpolation,
; then exponentiating on return, take log of grid points and function
; for interpolation of power-law behaviour.
; HISTORY:
; Written: Frank Varosi, HSTX @ NASA/GSFC, 1997.
; F.V. 1998, added /EXPONENTIAL_INTERP and /POWER_LAW_INTERP options.
; F.V. 1999, /EXP and /POWER options are remebered for repeat calls,
;   and scalar value is returned if result has just one element.
;-

function finterpol, fxold, xold, xnew, QUIET=quiet, INITIALIZE=init, $
          POWER_LAW_INTERP=ipow, $
          EXPONENTIAL_INTERP=iexp
   common finterpol, Lw, Lw1, xfrac
   common finterpol2, logflag, nextrap

  Nxo2 = N_elements( xold )-2

  if (Nxo2 ge 0) and (N_elements( xnew ) gt 0) then begin

    Lw = InterLeave( xold, xnew )
    wex = where( (Lw LT 0) or (Lw GT Nxo2), nextrap )

    if nextrap GT 0 then begin
      if Lw(wex(nextrap-1)) eq (Nxo2+1) then begin
        w = wex(nextrap-1)
        if xnew(w) eq xold(Lw(w)) then nextrap=nextrap-1
         endif
      Lw(wex) = ( Lw(wex) > 0 ) < Nxo2
       endif

    Lw1 = Lw+1
    if keyword_set( ipow ) or keyword_set( iexp ) then logflag=1 $
                else logflag=0
    if keyword_set( ipow ) then begin
      z = alog( xold > 1e-37 )
      xfrac = ( alog(xnew>1e-37) - z(Lw) )/( z(Lw1) - z(Lw) )
     endif else xfrac = ( xnew - xold(Lw) )/( xold(Lw1) - xold(Lw) )

   endif else if (Nxo2 eq -1) and (N_elements( xnew ) gt 0) then begin
    message,"need more than 1 point in the Xold grid",/INFO
    return,0
    endif

  if keyword_set( init ) then return,1

  if (N_elements( fxold ) LE 0) or (N_elements( xfrac ) LE 0) then begin
    if N_elements( fxold ) LE 0 then $
      message,"missing parameter fxold",/INFO
    if N_elements( xfrac ) LE 0 then $
      message,"missing parameters xold or xnew",/INFO
    print,"syntax:
    print," first call: fxnew = finterpol( fxold, xold, xnew )"
    print," repeat calls: fxnew = finterpol( fxold )
    return,0
    endif

  if (nextrap GT 0) and NOT keyword_set( quiet ) then  $
    message,"extrapolating at " + strtrim( nextrap, 2 ) + $
      " out of the " + strtrim( N_elements( xfrac ), 2 ) $
              + " new points",/INFO

  if keyword_set( logflag ) then begin

    z = alog( fxold > 1e-37 )
    result = exp( ( z(Lw1) - z(Lw) ) * xfrac + z(Lw) )

   endif else result = ( fxold(Lw1) - fxold(Lw) ) * xfrac + fxold(Lw)

  if N_elements( result ) eq 1 then return, scalar( result ) $
          else return, result
end
;--------------------------------------------------------------------------------


pro match, a, b, suba, subb, COUNT = count, SORT = sort
;+
; NAME:
;       MATCH
; PURPOSE:
;       Routine to match values in two vectors.
;
; CALLING SEQUENCE:
;       match, a, b, suba, subb, [ COUNT =, /SORT ]
;
; INPUTS:
;       a,b - two vectors to match elements, numeric or string data types
;
; OUTPUTS:
;       suba - subscripts of elements in vector a with a match
;               in vector b
;       subb - subscripts of the positions of the elements in
;               vector b with matchs in vector a.
;
;       suba and subb are ordered such that a[suba] equals b[subb]
;
; OPTIONAL INPUT KEYWORD:
;       /SORT - By default, MATCH uses two different algorithm: (1) the 
;               /REVERSE_INDICES keyword to HISTOGRAM is used for integer data,
;               while a sorting algorithm is used for non-integer data.   The
;               histogram algorithm is usually faster, except when the input
;               vectors are sparse and contain very large numbers, possibly
;               causing memory problems.   Use the /SORT keyword to always use
;               the sort algorithm.
;               
; OPTIONAL KEYWORD OUTPUT:
;       COUNT - set to the number of matches, integer scalar
;
; SIDE EFFECTS:
;       The obsolete system variable !ERR is set to the number of matches;
;       however, the use !ERR is deprecated in favor of the COUNT keyword 
;
; RESTRICTIONS:
;       The vectors a and b should not have duplicate values within them.
;       You can use rem_dup function to remove duplicate values
;       in a vector
;
; EXAMPLE:
;       If a = [3,5,7,9,11]   & b = [5,6,7,8,9,10]
;       then 
;               IDL> match, a, b, suba, subb, COUNT = count
;
;       will give suba = [1,2,3], subb = [0,2,4],  COUNT = 3
;       and       a[suba] = b[subb] = [5,7,9]
;
; 
; METHOD:
;       For non-integer data types, the two input vectors are combined and
;       sorted and the consecutive equal elements are identified.   For integer
;       data types, the /REVERSE_INDICES keyword to HISTOGRAM of each array
;       is used to identify where the two arrays have elements in common.   
; HISTORY:
;       D. Lindler  Mar. 1986.
;       Fixed "indgen" call for very large arrays   W. Landsman  Sep 1991
;       Added COUNT keyword    W. Landsman   Sep. 1992
;       Fixed case where single element array supplied   W. Landsman Aug 95
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use a HISTOGRAM algorithm for integer vector inputs for improved 
;             performance                W. Landsman         March 2000
;       Work again for strings           W. Landsman         April 2000
;       Use size(/type)                  W. Landsman         December 2002
;       Work for scalar integer input    W. Landsman         June 2003
;       Assume since V5.4, use COMPLEMENT to WHERE() W. Landsman Apr 2006
;-
;-------------------------------------------------------------------------
 On_error,2
 compile_opt idl2

 if N_params() LT 3 then begin
     print,'Syntax - match, a, b, suba, subb, [ COUNT = ]'
     print,'    a,b -- input vectors for which to match elements'
     print,'    suba,subb -- output subscript vectors of matched elements'
     return
 endif

 da = size(a,/type) & db =size(b,/type)
 if keyword_set(sort) then hist = 0b else $
 hist = (( da LE 3 ) or (da GE 12)) and  ((db LE 3) or (db GE 12 )) 

 if not hist then begin           ;Non-integer calculation
 
 na = N_elements(a)              ;number of elements in a
 nb = N_elements(b)             ;number of elements in b

; Check for a single element array

 if (na EQ 1) or (nb EQ 1) then begin
        if (nb GT 1) then begin
                subb = where(b EQ a[0], nw)
                if (nw GT 0) then suba = replicate(0,nw) else suba = [-1]
        endif else begin
                suba = where(a EQ b[0], nw)
                if (nw GT 0) then subb = replicate(0,nw) else subb = [-1]
        endelse
        count = nw
        return
 endif
        
 c = [ a, b ]                   ;combined list of a and b
 ind = [ lindgen(na), lindgen(nb) ]       ;combined list of indices
 vec = [ bytarr(na), replicate(1b,nb) ]  ;flag of which vector in  combined 
                                         ;list   0 - a   1 - b

; sort combined list

 sub = sort(c)
 c = c[sub]
 ind = ind[sub]
 vec = vec[sub]

; find duplicates in sorted combined list

 n = na + nb                            ;total elements in c
 firstdup = where( (c EQ shift(c,-1)) and (vec NE shift(vec,-1)), Count )

 if Count EQ 0 then begin               ;any found?
        suba = lonarr(1)-1
        subb = lonarr(1)-1
        return
 end
 
 dup = lonarr( Count*2 )                     ;both duplicate values
 even = lindgen( N_elements(firstdup))*2     ;Changed to LINDGEN 6-Sep-1991
 dup[even] = firstdup
 dup[even+1] = firstdup+1
 ind = ind[dup]                         ;indices of duplicates
 vec = vec[dup]                         ;vector id of duplicates
 subb = ind[ where( vec, complement = vzero) ]             ;b subscripts
 suba = ind[ vzero] 
  
 endif else begin             ;Integer calculation using histogram.

 minab = min(a, MAX=maxa) > min(b, MAX=maxb) ;Only need intersection of ranges
 maxab = maxa < maxb

;If either set is empty, or their ranges don't intersect: 
;  result = NULL (which is denoted by integer = -1)
  !ERR = -1
  suba = -1
  subb = -1
  COUNT = 0L
 if (maxab lt minab) or (maxab lt 0) then return
 
 ha = histogram([a], MIN=minab, MAX=maxab, reverse_indices=reva)
 hb = histogram([b], MIN=minab, MAX=maxab, reverse_indices=revb)
 
 r = where((ha ne 0) and (hb ne 0), count)
 if count gt 0 then begin
  suba = reva[reva[r]]
  subb = revb[revb[r]]
 endif 
 endelse 

 return
 
 end
;------------------------------------------------------------------------------
function notind, index, nelts
;+
; NAME:
;   NOTIND
; PURPOSE:
;   To generate the complementary indices to an array of index values
;   up to a given maximum.
;
; CALLING SEQUENCE:
;   not = NOTIND(index, n_elements)
;
; INPUTS:
;   INDEX -- Array of indices to generate the complement for.
;   N_ELEMENTS -- Number of elements out of which the inidices are coming.
; KEYWORD PARAMETERS:
;   NONE
;
; OUTPUTS:
;   NOT -- The complementary indices.
;
; MODIFICATION HISTORY:
;       Documented.
;       Wed Nov 21 12:18:52 2001, Erik Rosolowsky <eros@cosmic>
;-
x = bytarr(nelts)+1b
x[index] = 0
  return, where(x)
end

;--------------------------------------------------------------------------------
	FUNCTION TRIM, NUMBER, FORMAT, FLAG, QUIET=QUIET
;+
; Name        : 
;	TRIM()
; Purpose     : 
;	Converts numbers to strings, without trailing zeros.
; Explanation : 
;	Converts numbers into a string representation, and trims off leading
;	and/or trailing blanks.  Differs from STRTRIM in that trailing zeros
;	after the period are also trimmed off, unless NUMBER is already a
;	string, or an explicit format is passed.
; Use         : 
;	Result = TRIM( NUMBER  [, FORMAT ]  [, FLAG ] )
; Inputs      : 
;	NUMBER	= Variable or constant.  May be of any ordinary including
;		  string.  However, structures are not allowed.
; Opt. Inputs : 
;	FORMAT	- Format specification for STRING function.  Must be a string
;		  variable, start with the "(" character, end with the ")"
;		  character, and be a valid FORTRAN format specification.  If
;		  NUMBER is complex, then FORMAT will be applied separately to
;		  the real and imaginary parts.
;
;	FLAG	- Flag passed to STRTRIM to control the type of trimming:
;
;			FLAG = 0	Trim trailing blanks.
;			FLAG = 1	Trim leading blanks.
;			FLAG = 2	Trim both leading and trailing blanks.
;
;		  The default value is 2.  If NUMBER is complex, then FORMAT
;		  will be applied separately to the real and imaginary parts.
;
; Outputs     : 
;	Function returns as a string variable representing the value NUMBER.
; Opt. Outputs: 
;	None.
; Keywords    : 
;	None.
; Calls       : 
;	None.
; Common      : 
;	None.
; Restrictions: 
;	NUMBER must not be a structure.
;	FORMAT must be a valid format specification, and must not be passed
;		if NUMBER is of type string.
;	FLAG must not be of string type, or an array.
; Side effects: 
;	None.
; Category    : 
;	Utilities, Strings.
; Prev. Hist. : 
;	William Thompson	Applied Research Corporation
;	May, 1987		8201 Corporate Drive
;				Landover, MD  20785
;
;	William Thompson, Feb. 1992, added support for complex numbers, and
;				     fixed Unix problem with lowercase "e".
; Written     : 
;	William Thompson, GSFC, May 1987.
; Modified    : 
;	Version 1, William Thompson, GSFC, 9 April 1993.
;		Incorporated into CDS library.
;	Version 2, Zarro (SAC/GSFC), 3-Jun-98
;		Added check for undefined input
;       Version 3, Zarro (SM&A/GSFC), 1-Dec-99
;               Returned invalid input as blank string
;               to avoid downstream problems.
;       Version 4, Zarro (SM&A/GSFC), 4-Jan-00
;               Added /QUIET
;       Version 5, Zarro (SM&A/GSFC), 20-Jan-00
;               Vectorized
;	Version 6, 24-Jan-2000, William Thompson, GSFC
;		Fixed bug introduced in version 5.
;	Version 7, 14-Mar-2000, Zarro (SM&A/GSFC)
;		Moved check for unsupported type ahead of recursion
;-

	ON_ERROR,2
        loud=1-keyword_set(quiet)

;
;  Check for undefined input
;
        IF N_ELEMENTS(NUMBER) EQ 0 THEN BEGIN
         if loud then message,'Undefined input argument',/cont
         return,''
        ENDIF

;
;  Check the type of the variable NUMBER.
;
	S = SIZE(NUMBER)
	TYPE = S[S[0] + 1]
        IF TYPE GE 8 THEN BEGIN
         if loud then message,'Unsupported input argument',/cont
         return,''
        ENDIF

;
;  Check for vector input.
;
	IF S[0] GT 0 THEN BEGIN
	    NP = N_ELEMENTS(NUMBER)
	    IF NP GT 1 THEN OUT = MAKE_ARRAY(/STRING,DIMENSION=S[1:S[0]]) $
		    ELSE OUT = 'String'
	    CASE N_PARAMS() OF
		1: FOR I=0,NP-1 DO OUT[I] = TRIM(NUMBER[I],QUIET=QUIET)
		2: FOR I=0,NP-1 DO OUT[I] = TRIM(NUMBER[I],QUIET=QUIET,FORMAT)
		3: FOR I=0,NP-1 DO OUT[I] = TRIM(NUMBER[I],QUIET=QUIET,FORMAT,$
			FLAG)
	    ENDCASE
	    RETURN,OUT
	ENDIF
;
;  If NUMBER is complex, then process the real and imaginary parts separately.
;
	IF TYPE EQ 6 THEN BEGIN
		RNUMBER = FLOAT(NUMBER)
		INUMBER = IMAGINARY(NUMBER)
		CASE N_PARAMS() OF
			1:  BEGIN
				RNUMBER = TRIM(RNUMBER)
				INUMBER = TRIM(INUMBER)
				END
			2:  BEGIN
				RNUMBER = TRIM(RNUMBER,FORMAT)
				INUMBER = TRIM(INUMBER,FORMAT)
				END
			3:  BEGIN
				RNUMBER = TRIM(RNUMBER,FORMAT,FLAG)
				INUMBER = TRIM(INUMBER,FORMAT,FLAG)
				END
		ENDCASE
		RETURN, '(' + RNUMBER + ',' + INUMBER + ')'
	ENDIF
;
;  If only NUMBER was passed, then return the desired result.
;
	IF N_PARAMS(0) EQ 1 THEN BEGIN
		IF TYPE EQ 7 THEN BEGIN
			TRM = STRTRIM(NUMBER,2)
			GOTO,RETURN
		END ELSE BEGIN
			IF NUMBER EQ 0 THEN BEGIN
				TRM = '0'
				GOTO,RETURN
			END ELSE BEGIN
				TRM = STRTRIM( STRING(NUMBER), 2 )
				GOTO,REMOVE
			ENDELSE
		ENDELSE
	ENDIF
;
;  Check the type of the variable FORMAT.
;
	S = SIZE(FORMAT)
	TYPE_FORMAT = S[S[0] + 1]
;
;  If only two parameters were passed, then decide whether FORMAT or FLAG was
;  passed, and return the desired result. 
;
	IF N_PARAMS(0) EQ 2 THEN BEGIN
		IF TYPE_FORMAT EQ 7 THEN BEGIN
			TRM = STRTRIM( STRING(NUMBER,FORMAT), 2 )
			GOTO,RETURN
		END ELSE BEGIN
			FLAG = FORMAT
			IF TYPE EQ 7 THEN BEGIN
				TRM = STRTRIM( NUMBER, FLAG )
				GOTO,RETURN
			END ELSE BEGIN
				IF NUMBER EQ 0 THEN BEGIN
					TRM = '0'
					GOTO,RETURN
				END ELSE BEGIN
					TRM = STRTRIM( STRING(NUMBER), FLAG )
					GOTO,REMOVE
				ENDELSE
			ENDELSE
		ENDELSE
	ENDIF
;
;  All parameters were passed.  Act accordingly.
;
	TRM = STRTRIM( STRING(NUMBER,FORMAT), FLAG )
	GOTO,RETURN
;
;  Remove any trailing zeros.  First, check to make sure that the string 
;  contains a period.
;
REMOVE:
	TRM = STRUPCASE(TRM)
	IF STRPOS(TRM,'.') EQ -1 THEN GOTO,RETURN
;
;  Find and remove any exponential.
;
	LEN = STRLEN(TRM)
	EXP_POS = STRPOS(TRM,'E')
	IF EXP_POS EQ -1 THEN EXP = '' ELSE BEGIN
		EXP = STRMID(TRM,EXP_POS,LEN)
		TRM = STRMID(TRM,0,EXP_POS)
		LEN = STRLEN(TRM)
	ENDELSE
;
;  Keep removing trailing zeros until done.
;
	WHILE STRMID(TRM,LEN-1,1) EQ '0' DO BEGIN
		TRM = STRMID(TRM,0,LEN-1)
		LEN = LEN - 1
	ENDWHILE
;
;  If the last character is a period, remove it as well.
;
	IF STRMID(TRM,LEN-1,1) EQ '.' THEN BEGIN
		TRM = STRMID(TRM,0,LEN-1)
		LEN = LEN - 1
	ENDIF
;
;  Restore the exponential.
;
	TRM = TRM + EXP
;
;  Return the trimmed string TRM.
;
RETURN:
	RETURN,TRM
	END

;---------------------------------------------------------------------------
;+
;$Id: uniq_nosort.pro,v 1.1 2007/01/12 22:41:51 nathan Exp $
;
; Project     : SOHO - LASCO, STEREO - SECCHI 
;
; Name        : UNIQ_NOSORT
;
; Purpose     : Return the subscripts of the unique elements in an array.
;               Does not require array to be sorted (as in UNIQ).
;
; Use         : UNIQ_NOSORT(Array)
;
; Inputs      : Array  The array to be scanned.  
;
; Opt. Inputs : None
;
; Outputs     : An array of indicies into ARRAY is returned.  The expression:
;               ARRAY(UNIQ_NOSORT(ARRAY))
;               will be a copy of the sorted Array with duplicate elements removed.
;
; Opt. Outputs: None
;
; Keywords    : None
;
; Prev. Hist. : Adapted from SOHO/LASCO planning tool.
;
; Written     : Scott Paswaters, NRL, Dec 1996.
;
; @(#)uniq_nosort.pro   1.1 05/14/97 :NRL Solar Physics
;               
; Modification History:
;
; $Log: uniq_nosort.pro,v $
; Revision 1.1  2007/01/12 22:41:51  nathan
; moved from nrl_lib/dev/database
;
; Revision 1.1  2005/10/26 14:41:30  esfand
; Secchi database routines.
;
;--------------------------------------------------------------------------
;
FUNCTION UNIQ_NOSORT, a

   len = N_ELEMENTS(a)
   used = BYTARR(len)
   new = 0
   ind = WHERE(a EQ a(0))
   used(ind) = 1
   FOR i=1L, len-1 DO BEGIN
      IF (used(i) EQ 0) THEN BEGIN
         new = [new, i]
         ind = WHERE(a EQ a(i))
         used(ind) = 1
      ENDIF
   ENDFOR

   RETURN, new

END
;---------------------------------------------------;
;---------------SUPPORTING FUNCTIONS----------------;
;-----------------------END-------------------------;


;---------------------------------------------------;
;-----------------SVD_REG_MAIN----------------------;
;---------------------------------------------------;
FUNCTION SVD_REG_MAIN,x_norma,bbT,b_norm,GF_I_T,ch_ch,C_reg,dEp_b,ipmax,dEe_b,iemax,ifilt_edges,particles

; this is the NEW UNDER DEVELOPMENT CASE

; This is the main function of the srem_svd_noa code        


; UNDER DEVELOPMENT!

; Inputs:

; x_norma: the normalisation value for the (binned) proton and electron differential
; flux multiplied with the width of the corresponding energy bin
; bbT: vector with SREM count-rates from all 15 counters
; b_norm: The normalisation values of the count-rates
; GF_I_T: The response function re-binned for the selected proton and electron bins
; dEp_b,dEe_b; Energy widh of proton/electron energy bins
; ipmax,iemax; maximum electron/proton bin to be considered for the unfolding
; ch_ch; Chosen Channels (counters)
; C_reg: Regularisation matrix

; Output: Proton and electron omni-directional differential fluxes multiplied with the width of the corresponding energy bin

; I. Sandberg ISARS/NOA sandberg@noa.gr

FORWARD_FUNCTION REGSVD,REGSVD_ERR

    nx=iemax+ipmax+2
    nx_p=ipmax+1
    nx_e=iemax+1

     bb=bbT[ch_ch]         ; array of count-rates for selected counters
     Bcov= DIAG_MATRIX(bb) ; covariance matrix of count-rates

; Derivation of un-regularised solution

; Test lines: t - the x00 value should be the same with the one that follows
;    svdc,GF_I_T[*,ch_ch],Sini,U,V,/double  ; SVD for unregularised system
;    dummy=where(Sini gt 0,index)
;    w00=svsol(U,Sini,V,bb)


; We calculate the inverse of C_reg using the SVD components of C_reg
    svdc,C_reg,S_reg,U_reg,V_reg,/double ;

;    S_reg_inv=IMSL_INV(diag_matrix(S_reg),/DOUBLE)

    S_reg_inv=diag_matrix(1.0D/S_reg)


    C_inv=V_reg##S_reg_inv##transpose(U_reg);


    svdc,GF_I_T[*,ch_ch]##C_inv,S,U,V,/double ;
    d0=transpose(U)##bb

    w0=REGSVD(C_inv,V,d0,s,0.d0) ; unregularised soluition

  ; tau scanning ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  


          i=-1 
          repeat begin; 
          i=i+1


          w=REGSVD(C_inv,V,d0,s,s[i])

;          FPIO=w[0:nx_p-1-ifilt_edges] 
;          FEIO=w[nx_p:nx-1-ifilt_edges]


          FPIO=w[0:nx_p-1-ifilt_edges]*x_norma[0:ipmax-ifilt_edges] 
          FEIO=w[nx_p:nx-1-ifilt_edges]*x_norma[nx_p:nx-1-ifilt_edges]



          FPDO_diff=ts_diff(w[0:nx_p-1-ifilt_edges]*x_norma[0:ipmax-ifilt_edges]/(dEp_b[0:ipmax-ifilt_edges]),1)
          FEDO_diff=ts_diff(w[nx_p:nx-1-ifilt_edges]*x_norma[nx_p:nx-1-ifilt_edges]/(dEe_b[0:iemax-ifilt_edges]),1)


          dummy=where(FPIO lt 0.d0,ixxp)
          dummy=where(FEIO lt 0.d0,ixxe)

          dummy=where(FPDO_diff lt 0.d0,ixxpd)
          dummy=where(FEDO_diff lt 0.d0,ixxed)

;          condition_p=ixxp +ixxpd
;          condition_e=ixxe ;+ixxed


          endrep until ixxp*ixxe ne 0 or i eq n_elements(s)-2

         tau=s[i] ; The largest tau=s[i] that gives negative fluxes

         ivalid_0=0


if particles eq 'protons' then begin

         repeat begin

         tau=tau*1.1d0
        
          w=REGSVD(C_inv,V,d0,s,tau)


          FPIO=w[0:nx_p-1-ifilt_edges] 
          FEIO=w[nx_p:nx-1-ifilt_edges]

          FPIO_diff=ts_diff(w[0:nx_p-1-ifilt_edges]*x_norma[0:ipmax-ifilt_edges],1)
          FEIO_diff=ts_diff(w[nx_p:nx-1-ifilt_edges]*x_norma[nx_p:nx-1-ifilt_edges],1)

          FPDO_diff=ts_diff(w[0:nx_p-1-ifilt_edges]*x_norma[0:ipmax-ifilt_edges]/(dEp_b[0:ipmax-ifilt_edges]),1)
          FEDO_diff=ts_diff(w[nx_p:nx-1-ifilt_edges]*x_norma[nx_p:nx-1-ifilt_edges]/(dEe_b[0:iemax-ifilt_edges]),1)


          dummy=where([FPIO,FEIO] lt 0.d0,ixx)

          dummy=where(FPIO_diff lt 0.d0,ixxp)
          dummy=where(FEIO_diff lt 0.d0,ixxe)

          dummy=where(FPDO_diff lt 0.d0,ixxpd)
          dummy=where(FEDO_diff lt 0.d0,ixxed)

endrep until ixx+ixxpd+ixxed eq 0 

endif



if particles eq 'electrons' then begin

         repeat begin
         tau=tau*1.1d0
      
         w=REGSVD(C_inv,V,d0,s,tau)


          FPIO=w[0:nx_p-1-ifilt_edges] 
          FEIO=w[nx_p:nx-1-ifilt_edges]

          dummy=where([FPIO,FEIO] lt 0.d0,ixx)

          FPIO_diff=ts_diff(w[0:nx_p-1-ifilt_edges]*x_norma[0:ipmax-ifilt_edges],1)
          FEIO_diff=ts_diff(w[nx_p:nx-1-ifilt_edges]*x_norma[nx_p:nx-1-ifilt_edges],1)

          FPDO_diff=ts_diff(w[0:nx_p-1-ifilt_edges]*x_norma[0:ipmax-ifilt_edges]/(dEp_b[0:ipmax-ifilt_edges]),1)
          FEDO_diff=ts_diff(w[nx_p:nx-1-ifilt_edges]*x_norma[nx_p:nx-1-ifilt_edges]/(dEe_b[0:iemax-ifilt_edges]),1)


          dummy=where(FPIO_diff lt 0.d0,ixxp)
          dummy=where(FEIO_diff lt 0.d0,ixxe)

          dummy=where(FPDO_diff lt 0.d0,ixxpd)
          dummy=where(FEDO_diff lt 0.d0,ixxed)



endrep until ixx+ixxed+ixxpd eq 0 

endif

tauf=tau

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; derivation of regularised (normalised) solution
w=REGSVD(C_inv,V,d0,s,tauf) ;  normalised solution
Wcov=REGSVD_ERR(C_inv,V,d0,s,tauf) ; covariance matrix

Xcov=DIAG_MATRIX(x_norma)##Wcov##DIAG_MATRIX(x_norma) ; covariance matrix

x=x_norma*w                       ; solution (not normalised)
x_error=diag_matrix(sqrt(Xcov))   ; flux errors

iflx_p=x[0:nx_p-1-ifilt_edges]
iflx_e=x[nx_p:nx-1-ifilt_edges]     
iflx_p_err=x_error[0:nx_p-1-ifilt_edges]
iflx_e_err=x_error[nx_p:nx-1-ifilt_edges]     

SVD_reg_main_struct=create_struct('iflx_p',iflx_p,'iflx_p_err',iflx_p_err,'iflx_e',iflx_e,'iflx_e_err',iflx_e_err)
RETURN,SVD_reg_main_struct

END
;---------------------------------------------------;
;-----------------REG_MATRIX------------------;
;---------------------------------------------------;
FUNCTION REG_MATRIX,nx,nx_p
; construction of regularisation matrix
; I. Sandbeerg


 ksi=10.d0^(-2.d0) ; small diagonal component needed FOR the inversion
 C_reg=dblarr(nx,nx)
 FOR iy=0,nx-1 DO BEGIN & FOR ix=0,nx-1 DO BEGIN 
 IF iy eq ix THEN   C_reg[ix,iy]= -2.d0 + ksi 
 IF abs(iy-ix) eq 1 THEN C_reg[ix,iy]=1.d0 & ENDFOR & ENDFOR

  
  C_reg[0]=-1.d0 + ksi  &  C_reg[nx*nx-1]=-1.d0 + ksi  &   C_reg[nx_p-1,nx_p-1] = -1.d0+ksi  ; first and last points (r)

 ; last 'row of proton spectra

 IF nx ne nx_p THEN BEGIN
  C_reg[nx_p,nx_p-1] = 0.d0 & C_reg[nx_p,nx_p]  = -1.d0+ksi        &    C_reg[nx_p-1,nx_p] = 0.d0       ; first 'line' of electron spectra
  C_reg[nx_p:nx-1,*]=C_reg[nx_p:nx-1,*]
 ENDIF

  RETURN, C_reg
END

;---------------------------------------------------;
;--------------SREM_SVD_CALC_FLUX-------------;
;---------------------------------------------------;
FUNCTION SREM_SVD_CALC_FLUX, bb,GF_I_T_B,Ep_b,dEp_b,Ee_b,dEe_b,ch_ch,i_filt_edges,particles

;; This function calls the SVD_REG_MAIN function over different ranges of particles
; energy spectrum and constructs a final solution by choosing the maximum values of the derived spectra

;; It rerutns the Omnidirectional differential spectra: particles/cm^2/sec/MeV ###

FORWARD_FUNCTION SVD_REG_MAIN

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; INPUT:
;  bb:           vector of all SREM count-rates to be unfolded
;  GF_I_T_B:     response function in the rebinned energies
;  Ep_b, Ee_b:   vector with the log-centers of the proton and electron energy bins
;  dEp_b,dEe_b:  vector with the widths of proton and electron energy bins
;  ch_ch:        vector with the chosen channels to be considered in the unfolding 
; particles:     string variable that defines which spectra to be unfolded and provide output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  OUTPUT:
; structure fluxes
; that contains the "particles" differential flux in units [#/cm^2/sec/MeV], 
; the corresponding uncertaities and the associated reconstructed particle count-rates.
; if particles='protons' the output is

;fluxes =create_struct('fp_final',fp_final,'fp_final_err',fp_final_err,'bb_svd_max_p',bb_svd_max_p)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FORWARD_FUNCTION UNFOLD_FUNCTIONS, SVD_REG_MAIN, VANISH_EDGES, EXTRAP_NANS_POWER, FLX2CNTS, NAN2ZERO

nx_pp=n_elements(Ep_b)
nx_ee=n_elements(Ee_b)
  
      ipmin=1+i_filt_edges
      iemin=1+i_filt_edges
 
      ratio_counts_p=DBLARR(nx_pp-ipmin)
      ratio_counts_e=DBLARR(nx_ee-iemin)
      
      khi2_p=DBLARR(nx_pp-ipmin)
      khi2_e=DBLARR(nx_ee-iemin)
      
      
      ifp_f_all=DBLARR(nx_pp-ipmin,nx_pp)+!VALUES.F_NAN
      ife_f_all=DBLARR(nx_ee-iemin,nx_ee)+!VALUES.F_NAN
      ifp_f_all_err=DBLARR(nx_pp-ipmin,nx_pp)+!VALUES.F_NAN
      ife_f_all_err=DBLARR(nx_ee-iemin,nx_ee)+!VALUES.F_NAN

      fe_f_all=DBLARR(nx_ee-iemin,nx_ee)+!VALUES.F_NAN
      fp_f_all=DBLARR(nx_pp-ipmin,nx_pp)+!VALUES.F_NAN      



if particles eq 'protons' then begin

      ; proton-do-loop ;

      iemax=nx_ee-1

      FOR ipmax= ipmin,nx_pp-1   DO BEGIN

         nx=iemax+ipmax+2
         nx_p=ipmax+1
         nx_e=iemax+1

         x_norma=DBLARR(nx)+1.d0 ; normalisation flux
         b_norm=DBLARR(15)+1
         b_norm=SQRT(ABS(bb))

         GF_I_T=[GF_I_T_B[0:ipmax,*],GF_I_T_B[nx_pp:nx_pp+iemax,*]]
         
         
         FOR i=0,nx-1 DO BEGIN &  GF_I_T[i,*]=GF_I_T[i,*]*x_norma[i] & ENDFOR
         FOR j=0,14 DO BEGIN & GF_I_T[*,j]=GF_I_T[*,j]/b_norm[j] & ENDFOR
       
  ;RESCALING AND ROTATIONS
         arr=TRANSPOSE([dEp_b[0:ipmax],dEe_b[0:iemax]]#(INTARR(nx)+1))
         C_reg=REG_MATRIX(nx,nx_p)
         C_reg=C_reg*arr   ;/max([dEp_b,dEe_b/max(dEp_b)])
         C_reg=C_reg/MAX(C_reg)
         bbT=bb[ch_ch]/b_norm[ch_ch]
         GF_I_Tn=GF_I_T
         
    
              
         FOR i=0,nx-1 DO BEGIN &  GF_I_Tn[i,*]=GF_I_T[i,*]*x_norma[i] & ENDFOR


         SVD_reg_main_struct=SVD_REG_MAIN(x_norma,bbT,b_norm,GF_I_Tn,ch_ch,C_reg,dEp_b,ipmax,dEe_b,iemax,i_filt_edges,particles)

                                          
         ifp_f=SVD_reg_main_struct.iflx_p
         ifp_f_err=SVD_reg_main_struct.iflx_p_err 
         
         npd=N_ELEMENTS(ifp_f)
         ifp_f=VANISH_EDGES(ifp_f,Ep_b[0:npd-1]) ; set each of flux solutions equlat to NaN in case of local minima
         ifp_f_all[ipmax-ipmin,0:npd-1]=ifp_f
         ifp_f_all_err[ipmax-ipmin,0:npd-1]=ifp_f_err
;         fp_f_all[ipmax-ipmin,0:npd-1]=ifp_f/dEp_b
ENDFOR ;npmax
      
;;;;;;;;; END OF FLUX DERIVATIONS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Synthetic solutions
ifp_final=(max(ifp_f_all,maxis,dimension=1,/nan))
ifp_final_err=ifp_f_all_err[maxis]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ifp_final=EXTRAP_NANS_POWER(ifp_final/dEp_b,Ep_b)*dEp_b

fp_final=ifp_final/dEp_b
fp_final_err=EXTRAP_NANS_POWER((ifp_final+ifp_final_err)/dEp_b,Ep_b)-fp_final

bb_svd_max_p=FLX2CNTS(NAN2ZERO(ifp_final,0.d0)/x_norma[0:nx_p-1],GF_I_T[0:nx_p-1,*])
bb_svd_max_p=bb_svd_max_p*b_norm

fluxes=create_struct('fp_final',fp_final,'fp_final_err',fp_final_err,'bb_svd_max_p',bb_svd_max_p)
; fluxes=create_struct('fp_final',fp_final,'fp_final_err',fp_final_err,'bb_svd_max_p',bb_svd_max_p,'fp_f_all',fp_f_all)

      ;END of proton loop
RETURN,fluxes      
ENDIF


if particles eq 'electrons' then begin
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; electron loop
      ipmax= nx_pp-1
      FOR iemax= iemin,nx_ee-1   DO BEGIN ; do loop for width of energy spectra to be unfolded
         nx=iemax+ipmax+2
         nx_p=ipmax+1
         nx_e=iemax+1
         x_norma=DBLARR(nx)+1.d0 ; normalisation flux
         b_norm=DBLARR(15)+1
         b_norm=SQRT(ABS(bb))
         GF_I_T=[GF_I_T_B[0:ipmax,*],GF_I_T_B[nx_pp:nx_pp+iemax,*]]
         FOR i=0,nx-1 DO BEGIN &  GF_I_T[i,*]=GF_I_T[i,*]*x_norma[i] & ENDFOR
         FOR j=0,14 DO BEGIN & GF_I_T[*,j]=GF_I_T[*,j]/b_norm[j] & ENDFOR
         
         ;RESCALING AND ROTATIONS
         arr=TRANSPOSE([dEp_b[0:ipmax],dEe_b[0:iemax]]#(INTARR(nx)+1))
         C_reg=REG_MATRIX(nx,nx_p)
         C_reg=C_reg*arr;/max([dEp_b,dEe_b/max(dEp_b)])
         C_reg=C_reg/MAX(C_reg)
         
         bbT=bb[ch_ch]/b_norm[ch_ch]
         GF_I_Tn=GF_I_T
         
         FOR i=0,nx-1 DO BEGIN &  GF_I_Tn[i,*]=GF_I_T[i,*]*x_norma[i] & ENDFOR
         

          SVD_reg_main_struct=SVD_REG_MAIN(x_norma,bbT,b_norm,GF_I_T,ch_ch,C_reg,dEp_b,ipmax,dEe_b,iemax,i_filt_edges,particles)
 

         ife_f=SVD_reg_main_struct.iflx_e
         ife_f_err=SVD_reg_main_struct.iflx_e_err
         ned=N_ELEMENTS(ife_f)
         ife_f=VANISH_EDGES(ife_f,Ee_b[0:ned-1])
         ife_f_all[iemax-iemin,0:ned-1]=ife_f
         ife_f_all_err[iemax-iemin,0:ned-1]=ife_f_err
;         fe_f_all[iemax-iemin,0:ned-1]=ife_f/dEe_b
     ENDFOR ;nemax

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; END OF FLUX DERIVATIONS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Synthetic solutions
ife_final=(max(ife_f_all,maxis,dimension=1,/nan))
ife_final_err=ife_f_all_err[maxis]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




; Fill the NaN with power-law interpolation
ife_final=EXTRAP_NANS_POWER(ife_final/dEe_b,Ee_b)*dEe_b

; differential fluxes
fe_final=ife_final/dEe_b

fe_final_err=EXTRAP_NANS_POWER((ife_final+ife_final_err)/dEe_b,Ee_b)-fe_final

bb_svd_max_e=FLX2CNTS(NAN2ZERO(ife_final,0.d0)/x_norma[nx_p:nx-1],GF_I_T[nx_p:nx-1,*])
bb_svd_max_e=bb_svd_max_e*b_norm

; fluxes=create_struct('fe_final',fe_final,'fe_final_err',fe_final_err,'bb_svd_max_e',bb_svd_max_e,'fe_f_all',fe_f_all)

fluxes=create_struct('fe_final',fe_final,'fe_final_err',fe_final_err,'bb_svd_max_e',bb_svd_max_e)

RETURN,fluxes     
ENDIF


END

;---------------------------------------------------;
;-------------------error---------------------;
;---------------------------------------------------;

FUNCTION error,bb,bb_rec
; the error- khi2 function used to evaluate results

error1=(bb-bb_rec)^2.d0/(bb^2.d0+bb_rec^2.d0)
error=total(error1,1,/nan,/double)
RETURN,error
END

;---------------------------------------------------;
;---------------------INDEX_FLUXES------------------;
;---------------------------------------------------;
FUNCTION INDEX_FLUXES, bb_p_rec,bb_e_rec,cnts

COMPILE_OPT STRICTARR

; Label for the unfolded fluxes
;   1: NOT RELIABLE ELECTRONS - STRONG PROTON FLUX
;  -1: NOT RELIABLE PROTONS   - STRONG ELECTRON FLUX
;   0: NOT EVIDENT DOMINATION OF PARTICLES
; I. Sandberg, ISARS/NOA

nt=n_elements(cnts)/15

index=intarr(nt)

d13=[0,1,2,11,12,13] ; counters to be used FOR the index

khi2_D13 =error(bb_e_rec[d13,*]+bb_p_rec[d13,*],cnts[d13,*])
khi2p_D13=error(bb_p_rec[d13,*],cnts[d13,*])
khi2e_D13=error(bb_e_rec[d13,*],cnts[d13,*])

protons   =  where(khi2p_D13 LT khi2_D13,inp)
electrons =  where(khi2e_D13 LT khi2_D13,ine)

IF inp GT 0 then index[protons]=1
IF ine GT 0 then index[electrons]=-1

; when both proton and electron indices are applicable we set index equal to zero
match,protons,electrons,suba,subb
IF suba[0] NE -1 then strange=protons[suba]
IF suba[0] NE -1 then index[strange]=0

return, index
END


;---------------------------------------------------;
;------------------REGSVD---------------------;
;---------------------------------------------------;

FUNCTION REGSVD,C_inv,V,d0,s,tau
; SVD regularized solution
;
             zz=d0*s/(s^2.d0+tau^2.d0)
;             d=d0*s^2.d0/(s^2.d0+tau^2.d0)
             w=C_inv##V##zz ;           
RETURN,w
END

;---------------------------------------------------;
;------------------REGSVD_ERR-----------------;
;---------------------------------------------------;

FUNCTION REGSVD_ERR,C_inv,V,d0,s,tau
; Covariance Matrix of SVD regularized solution
;
Zcov=DIAG_MATRIX(s^2.d0/((s^2.d0+tau^2.d0)^2.d0))
Wcov=C_inv##V##Zcov##(transpose(V))##(transpose(C_inv))
RETURN,Wcov
END

;---------------------------------------------------;
;------------------FLX2CNTS-------------------;
;---------------------------------------------------;

FUNCTION FLX2CNTS, w_t, GF
; calculate the 15 countrates for
; w_t: given integral fluxes within specified energy bins
; in accordance to GF. OBS: W=f(E)ΔΕ
;
;
nt=n_elements(w_t)/(n_elements(GF)/15)
bb_t=dblarr(15,nt)

FOR it=0,nt-1 DO BEGIN
bb_t[*,it]=GF##w_t[*,it]
ENDFOR
RETURN,bb_t
END

;---------------------------------------------------;
;------------------VANISH_EDGES---------------;
;---------------------------------------------------;


FUNCTION VANISH_EDGES,f,E
;
; If local minima exists, I set equal to NaN the values after
; I. Sandberg
;

FORWARD_FUNCTION extrema
f1=f
n=n_elements(f1)

IF n lt 4 THEN goto,RETURNH

maxima=extrema(f1,/MAX_ONLY)
minima=extrema(f1,/MIN_ONLY)

maxima=min(maxima)
minima=max(minima)

; I DO NOT KEEP THE MINIMA POINT
IF minima gt 1 THEN f1[minima-1:n-1]=!Values.F_NAN

; I DO NOT KEEP THE MAXIMA POINT
;IF maxima gt 1 THEN f[0:maxima-1]=!Values.F_NAN
RETURNH:
RETURN,f1
END

;---------------------------------------------------;
;------------------NAN2ZERO-------------------;
;---------------------------------------------------;

FUNCTION NAN2ZERO,f,value
; replace NANs with a fixed "error" value
; I. Sandberg

f1=f
n=n_elements(f1)
arr=WHERE(FINITE(f1, /NAN))  
IF arr[0] gt -1 THEN BEGIN
f1[arr]=value
ENDIF
RETURN,f1
END

;---------------------------------------------------;
;------------------value2NaN------------------------;
;---------------------------------------------------;

FUNCTION value2NaN,f,value
; replace zeros with NANs
; I. Sandberg

n=n_elements(f)
arr=WHERE(f eq value,arr)
IF arr[0] gt -1 THEN BEGIN
f[arr]=!Values.F_NAN
ENDIF
RETURN,f
END

;---------------------------------------------------;
;------------------zero2NaN-------------------------;
;---------------------------------------------------;

FUNCTION zero2NaN,f
; replace zeros with NANs
; I. Sandberg

n=n_elements(f)
arr=WHERE(f eq 0,arr)
IF arr[0] gt -1 THEN BEGIN
f[arr]=!Values.F_NAN
ENDIF
RETURN,f
END

;---------------------------------------------------;
;-------------EXTRAP_NANS_POWER---------------;
;---------------------------------------------------;

FUNCTION EXTRAP_NANS_POWER,f,E
; replace NaNs with power-law interpolated values
; I. Sandberg

FORWARD_FUNCTION FINTERPOL

f1=f
n=n_elements(f1)
arr=WHERE(FINITE(f, /NAN),n_nans)

IF n_nans gt 0 THEN BEGIN
arrv=notind(arr,n)
arr=dblarr(n_nans)+arr

   IF n_nans gt 1 THEN f1[arr]=finterpol(f1[arrv],E[arrv],E[arr],/power_law_interp,/quiet)

   IF n_nans eq 1 THEN BEGIN
   dummy=FINTERPOL(f1[arrv],E[arrv],[E[arr],E[arr]],/power_law_interp,/quiet) ; cheap tryck IF interpolated is a scalar
   f1[arr]=dummy[0]
   ENDIF

ENDIF
RETURN,f1
END

;---------------------------------------------------;
;-------LINEAR_INTERP_SUCCESIVE---------------;
;---------------------------------------------------;
FUNCTION LINEAR_INTERP_SUCCESIVE, y,x,xnew

; This function locates the neighboring points
; and performs a piece-wise linear interpolation

n=n_elements(xnew)
ny=n_elements(y)/n_elements(x)
;
bounds=dblarr(2,n) ; 
y2=dblarr(n,ny)
;
FOR i=0,n-1 DO BEGIN
bounds(0,i)=max(where(x le xnew[i]))
bounds(1,i)=min(where(x gt xnew[i]))
ENDFOR
;
;
FOR iy=0,ny-1 DO BEGIN ; do-loop in levels/counters

   FOR i=0,n-1 DO BEGIN ; do-loop in new energy values
   y2[i,iy]=interpol(y[iy,bounds[*,i]],x(bounds[*,i]),xnew[i]);  

ENDFOR
ENDFOR
;
RETURN,y2
;
END

;---------------------------------------------------;
;----------------A_I_q_Ebins------------------;
;---------------------------------------------------;
FUNCTION A_I_q_Ebins, nx_p, Ep_tot,E_p_bou,y

; FUNCTION that calculates the response at requested energy
; values through "proper" averaging into the energy bins

y_b=dblarr(15,nx_p)

FOR iich=0,14 DO BEGIN ;
      FOR i_nx_p=0,nx_p-1 DO BEGIN
      ind1=where(Ep_tot eq E_p_bou[i_nx_p])
      ind2=where(Ep_tot eq E_p_bou[i_nx_p+1])     
;PRINT,Ep_tot,E_p_bou[i_nx_p],E_p_bou[i_nx_p+1]
;PRINT,IND1,IND2
      ff=y[iich,ind1:ind2] &   aa=Ep_tot[ind1:ind2]
      y_b[iich,i_nx_p]=INT_TABULATED(aa , ff, /DOUBLE , /SORT)/(Ep_tot[ind2]-Ep_tot[ind1])
      IF  y_b[iich,i_nx_p] lt 0 THEN y_b[iich,i_nx_p]=0.d0    
      ENDFOR
ENDFOR
  RETURN,y_b

END

;---------------------------------------------------;
;----------------GF_BINS----------------------------;
;---------------------------------------------------;
FUNCTION GF_BINS, i_srem_unit,Ep_b,Ee_b,Ep_bin_bou,Ee_bin_bou

;HELP, i_srem_unit,Ep_b,Ee_b,Ep_bin_bou,Ee_bin_bou

;STOP
; This function calculates the response - GF - for requested energy bins and SREM unit.
; The following actions are taken:
; (1) Restore calibration GF and calibration energies for selected SREM unit
; (2) A linear fit of GF at the edges of the requested energy bins by
; using the values of the first neighboring points of calibration energies
; (3) Average and find GF within each requested bin for the calculation of GF_output

; written by I Sandberg


FORWARD_FUNCTION A_I_q_Ebins, LINEAR_INTERP_SUCCESIVE

; define outputs
nx_p=n_elements(Ep_b)
nx_e=n_elements(Ee_b)
GF_I_p_b=dblarr(nx_p,15) 
GF_I_e_b=dblarr(nx_e,15) 

; restore GF's

IF i_srem_unit eq 0 THEN restore,'/home/sandberg/data/SREM_DATA/SREM/SREM_GF_GIOVEB.idat'
IF i_srem_unit eq 1 THEN restore,'/home/sandberg/data/SREM_DATA/SREM/SREM_GF_PROBA1.idat'
IF i_srem_unit eq 2 THEN restore,'/home/sremdc/esa/data/SREM_DATA/SREM/SREM_GF_INTEGRAL.idat'
IF i_srem_unit eq 3 THEN restore,'/home/sandberg/data/SREM_DATA/SREM/SREM_GF_ROSETTA.idat'
IF i_srem_unit eq 4 THEN restore,'/home/sandberg/data/SREM_DATA/SREM/SREM_GF_HERSCHEL.idat'
IF i_srem_unit eq 5 THEN restore,'/home/sandberg/data/SREM_DATA/SREM/SREM_GF_PLANCK.idat'

; energy for each (calibrating) energy bin
Epbin=transpose((E1_p_SREM_resp+E2_p_SREM_resp)/2.d0) &
Eebin=transpose((E1_e_SREM_resp+E2_e_SREM_resp)/2.d0)

; rewrite bin boundaries into a single vector
E_p_bou=Ep_bin_bou(uniq_nosort(Ep_bin_bou)) ; bin_boundaries
E_e_bou=Ee_bin_bou(uniq_nosort(Ee_bin_bou)) ; bin_boundaries

; calculate the GF on the edges of energy bins
GFp2=LINEAR_INTERP_SUCCESIVE(GF_P,Epbin,E_p_bou)
GFe2=LINEAR_INTERP_SUCCESIVE(GF_E,Eebin,E_e_bou)

; merge the selected and the calibrating energies and GFs 
; and we sort them in terms of energy ascENDing order

Ep_tot=[Epbin,E_p_bou]
Gp_tot=[[(GF_P)],[transpose(GFp2)]]
sortp=sort(Ep_tot)
Gp_tot=Gp_tot[*,sortp]
Ep_tot=Ep_tot[sortp]
;
Ee_tot=[Eebin,E_e_bou]
Ge_tot=[[(GF_E)],[transpose(GFe2)]]
sorte=sort(Ee_tot)
Ge_tot=Ge_tot[*,sorte]
Ee_tot=Ee_tot[sorte]
;
; calculate the GF's FOR the bins
GF_I_p_b=A_I_q_Ebins(nx_p, Ep_tot,E_p_bou,Gp_tot)
GF_I_e_b=A_I_q_Ebins(nx_e, Ee_tot,E_e_bou,Ge_tot)
GF_I_T=transpose([[(GF_I_p_b)],[(GF_I_e_b)]]) ; Total response
RETURN, GF_I_T
END

;---------------------------------------------------;
;------------ENERGY_FLUX_LEVELS---------------;
;---------------------------------------------------;

FUNCTION ENERGY_FLUX_LEVELS, Emins, Emaxs, nx_q, logspace=logspace,linspace=linspace,SEPEM=SEPEM
; For given energy range and number of bins calculate 
; (1) Bin Boundaries (2) Bin centers (3) Bin width
; I. Sandberg
; Options: logspace, linspace, SEPEM



IF keyword_set(SEPEM) THEN BEGIN
Eqmin=5.d0
Eqmax=200.d0
nx_q=10
dEq=(alog10(Eqmax)-alog10(Eqmin))/(2.d0*nx_q)             ; dEp of each half bin

nx_q=14
dEq_b=dblarr(nx_q) &  Eq_b=dblarr(nx_q) & Eq_bin_bou=dblarr(nx_q,2)
Eq15=10.d0^(alog10(Eqmin)+(findgen(2.d0*nx_q+1))*dEq)       ; E's that define the bin-boundaries and the bin centers

FOR iq=0,nx_q-1 DO BEGIN
Eq_bin_bou[iq,*]=[Eq15[2*(iq+1)-2],Eq15[2*(iq+1)]] ; BIN BOUNDARIES
dEq_b[iq]=Eq15[2*(iq+1)]-Eq15[2*(iq+1)-2]  &   Eq_b[iq]=Eq15[2*(iq+1)-1] ;
ENDFOR 

Eq_bin_bou=Eq_bin_bou[1:13,*]
Eq_bin_bou(0)=9.d0
Eq_bin_bou(n_elements(Eq_bin_bou)-1)=800.d0
dEq_b=dEq_b[1:13]
Eq_b=Eq_b[1:13]
dEq_b(0)=Eq_bin_bou[0,1]-Eq_bin_bou[0,0]
Eq_b(0)=Eq_bin_bou[0,0]+(Eq_bin_bou[0,1]-Eq_bin_bou[0,0])/2.d0
dEq_b(12)=Eq_bin_bou[12,1]-Eq_bin_bou[12,0]
Eq_b(12)=Eq_bin_bou[12,0]+(Eq_bin_bou[12,1]-Eq_bin_bou[12,0])/2.d0
ENDIF ; SEPEM



 IF keyword_set(logspace) THEN BEGIN
 Eqmin=min(Emins) & Eqmax=max(Emaxs)
 dEq=(alog10(Eqmax)-alog10(Eqmin))/(2.d0*nx_q)             ; dEp of each half bin
 dEq_b=dblarr(nx_q) &  Eq_b=dblarr(nx_q) & Eq_bin_bou=dblarr(nx_q,2)
 Eq15=10.d0^(alog10(Eqmin)+(findgen(2.d0*nx_q+1))*dEq)       ; E's that define the bin-boundaries and the bin centers
    FOR iq=0,nx_q-1 DO BEGIN
    Eq_bin_bou[iq,*]=[Eq15[2*(iq+1)-2],Eq15[2*(iq+1)]] ; BIN BOUNDARIES
    dEq_b[iq]=Eq15[2*(iq+1)]-Eq15[2*(iq+1)-2]  &   Eq_b[iq]=Eq15[2*(iq+1)-1] ;
    ENDFOR 
 ENDIF

  IF keyword_set(linspace) THEN BEGIN
  Eqmin=min(Emins) & Eqmax=max(Emaxs)
  dEq=((Eqmax)-(Eqmin))/(2.d0*nx_q)             ; dEp of each half bin
  dEq_b=dblarr(nx_q) &  Eq_b=dblarr(nx_q) & Eq_bin_bou=dblarr(nx_q,2)
  Eq15=((Eqmin)+(findgen(2.d0*nx_q+1))*dEq)       ; E's that define the bin-boundaries and the bin centers
    FOR iq=0,nx_q-1 DO BEGIN
    Eq_bin_bou[iq,*]=[Eq15[2*(iq+1)-2],Eq15[2*(iq+1)]] ; BIN BOUNDARIES
    dEq_b[iq]=Eq15[2*(iq+1)]-Eq15[2*(iq+1)-2]  &   Eq_b[iq]=Eq15[2*(iq+1)-1] ;
  ENDFOR 
 ENDIF
 
e_f_l=create_struct('Eq_bin_bou',Eq_bin_bou,'dEq_b',dEq_b,'Eq_b',Eq_b)
RETURN,e_f_l
;Eq_bin_bou=e_f_l.Eq_bin_bou
;dEq_b=e_f_l.dEq_b
;Eq_b=e_f_l.Eq_b
END


;------------------------------------------------------------------;
;-------EXTRACT SREM COUNTS FROM A PACC FILE-----------------------;
;------------------------------------------------------------------;
FUNCTION extract_counts, PACC_file

  id=CDF_OPEN(PACC_file)
  res=CDF_INQUIRE(id)
  CDF_CONTROL,id,GET_VAR_INFO=info,VARIABLE='COUNTRATE'
  for i=0,res.NZVARS-1 do begin
  res=CDF_VARINQ(id,i,/ZVARIABLE)
  CDF_CONTROL,id,VARIABLE=i,/ZVARIABLE,GET_VAR_INFO=V,SET_PADVALUE=0.0
  max_recs=V.MAXRECS
  endfor

  CDF_VARGET,id,'COUNTRATE', REC_COUNT=max_recs+1,  cnts,/zvariable ; 
  CDF_CLOSE,id

RETURN, cnts
END
;--------------------------------------------------------------------;
;------------Calculate SREM FLUXES from SREM COUNTS------------------;
;--------------------------------------------------------------------;
FUNCTION  COUNTS_TO_FLUXES, cnts, i_srem_unit

FORWARD_FUNCTION UNFOLD_FUNCTIONS, SREM_SVD_CALC_FLUX

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ch_ch=indgen(15) ;  list of counters to be used FOR the unfolding
i_filt_edges=2   ;  number of right edge flux levels to be ignored in the selection scheme of regularisation parameter
Epmax=400.d0     ;  maximum proton energy of spectra,   i.e. upper border of last energy proton bin
Eemax= 4.d0      ;  maximum electron energy of spectra, i.e. upper border of last energy electron bin
Epmin=11.d0      ;  minimum proton energy to be considered
;Eemin=0.55d0     ;  minimum electron energy to be considered
Eemin=0.65d0     ;  minimum electron energy to be considered

;
nx_pp_protons=15;  (>10)    number of bins covering the proton energy range [Epmin,Epmax]
nx_ee_protons=10;  (>7)
;
nx_pp_electrons=10;     number of bins covering the electron energy range [Eemin,Eemax]
nx_ee_electrons=15; 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Proton unfolding parameters for proton energy loop
;
nxx_protons=nx_pp_protons+nx_ee_protons

; Logarithmic discretisation of energy ranges & calculation of Geometric Factors
; for energy binning tuned for proton unfolding

pfl_protons=ENERGY_FLUX_LEVELS(Epmin,Epmax,nx_pp_protons,/logspace)
dEp_b_protons=pfl_protons.dEq_b & Ep_b_protons=pfl_protons.Eq_b & Ep_bin_bou_protons=pfl_protons.Eq_bin_bou
efl_protons=ENERGY_FLUX_LEVELS(Eemin,Eemax,nx_ee_protons,/logspace)
dEe_b_protons=efl_protons.dEq_b & Ee_b_protons=efl_protons.Eq_b & Ee_bin_bou_protons=efl_protons.Eq_bin_bou
GF_I_T_protons=GF_BINS(i_srem_unit,Ep_b_protons,Ee_b_protons,Ep_bin_bou_protons,Ee_bin_bou_protons)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Electron unfolding parameters for electron energy loop

nxx_electrons=nx_pp_electrons+nx_ee_electrons

; Logarithmic discetisation of energy ranges & calculation of Geometric Factors
; for energy binning tuned for electron unfolding

pfl_electrons=ENERGY_FLUX_LEVELS(Epmin,Epmax,nx_pp_electrons,/logspace)
dEp_b_electrons=pfl_electrons.dEq_b & Ep_b_electrons=pfl_electrons.Eq_b & Ep_bin_bou_electrons=pfl_electrons.Eq_bin_bou
efl_electrons=ENERGY_FLUX_LEVELS(Eemin,Eemax,nx_ee_electrons,/logspace)
dEe_b_electrons=efl_electrons.dEq_b & Ee_b_electrons=efl_electrons.Eq_b & Ee_bin_bou_electrons=efl_electrons.Eq_bin_bou
GF_I_T_electrons=GF_BINS(i_srem_unit,Ep_b_electrons,Ee_b_electrons,Ep_bin_bou_electrons,Ee_bin_bou_electrons)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The geometric functions of the derived flux levels - to be used in output files

GFP=GF_I_T_protons(0:nx_pp_protons-1,*)
GFE=GF_I_T_electrons(nx_pp_electrons:nx_pp_electrons+nx_ee_electrons-1,*)

nt=n_elements(cnts)/15
FPDO=dblarr(nx_pp_protons,nt)
FEDO=dblarr(nx_ee_electrons,nt)

FPDO_ERR=dblarr(nx_pp_protons,nt)
FEDO_ERR=dblarr(nx_ee_electrons,nt)

bb_p_rec=dblarr(15,nt)
bb_e_rec=dblarr(15,nt)

; Set a small number FOR zero measurements - needed FOR the rescaling
ICNTS0=where(cnts eq 0.d0)
IF ICNTS0[0] NE -1 THEN cnts[where(cnts eq 0.d0)]=0.0048875855d0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; D0-LOOP IN COUNT-RATE MEASUREMENTS (FOR OF EACH FILE)


FOR it=0,nt-1  DO BEGIN

bb=cnts[*,it] & IF total(finite(bb,/NAN)) GT 0 then goto,outputs

fluxes_protons=SREM_SVD_CALC_FLUX(bb,GF_I_T_protons,Ep_b_protons,dEp_b_protons,Ee_b_protons,dEe_b_protons,ch_ch,i_filt_edges,'protons') 


    FPDO[*,it]=fluxes_protons.fp_final
    FPDO_ERR[*,it]=fluxes_protons.fp_final_err
    bb_p_rec[*,it]=fluxes_protons.bb_svd_max_p

    fluxes_electrons=SREM_SVD_CALC_FLUX(bb,GF_I_T_electrons,Ep_b_electrons,dEp_b_electrons,Ee_b_electrons,dEe_b_electrons,ch_ch,i_filt_edges,'electrons')  


    FEDO[*,it]=fluxes_electrons.fe_final
    FEDO_ERR[*,it]=fluxes_electrons.fe_final_err
    bb_e_rec[*,it]=fluxes_electrons.bb_svd_max_e

outputs:
ENDFOR


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Derive index to label measurements

PARTICLE_INDEX=transpose(INDEX_FLUXES(bb_p_rec,bb_e_rec,cnts))


; convert fluxes to particles/MeV/cm^2/sec/str

   strrads=4.d0*!pi ; CONVERT EVERYTHING TO STEREORADs
   
   FPDO=FPDO/strrads
   FPDO_ERR=FPDO_ERR/strrads
   FEDO=FEDO/strrads
   FEDO_ERR=FEDO_ERR/strrads
   
   FPDO_ENERGY_bou=Ep_bin_bou_protons
   FEDO_ENERGY_bou=Ee_bin_bou_electrons

   FPDO_ENERGY=Ep_b_protons
   FEDO_ENERGY=Ee_b_electrons


PARTICLE_INDEX=NAN2ZERO(PARTICLE_INDEX,127)

FPDO    =transpose(NAN2ZERO(FPDO,fillval))
FEDO    =transpose(NAN2ZERO(FEDO,fillval))

help,fpdo,cnts,fedo

srem_fluxes=create_struct('FPDO',FPDO,'FEDO',FEDO,'FPDO_ENERGY',FPDO_ENERGY,'FEDO_ENERGY',FEDO_ENERGY,'FPDO_ENERGY_bou',FPDO_ENERGY_bou,'FEDO_ENERGY_bou',FEDO_ENERGY_bou,'PARTICLE_INDEX',PARTICLE_INDEX,'GFP',GFP,'GFE',GFE)

return,srem_fluxes


END

;--------------------------------------------------------------------;
;--!!!!!!---MAIN ALGORITHM-----LEVEL1-----TO-----LEVEL2----!!!!!-----;
;--------------------------------------------------------------------;
pro ESA_LEVEL1_to_LEVEL2
;
; Requires the functions:
; extract_counts.pro
; COUNTS_TO_FLUXES.pro
;-------------------------------------------------;
;---Input LEVEL 1: L1_cdf_file_name (PACC file)---;
;-------------------------------------------------;

PACC_file='/home/sremdc/esa/input/IREM_PACC_20050117_ESA_L1_P01.cdf'
I=0
IF (((I = STRPOS(PACC_file, 'SREMGIOVEB_'))) NE -1) THEN BEGIN
sat='GIOVEB'
i_srem_unit=0
ENDIF

IF (((I = STRPOS(PACC_file, 'SREMPROBA1_'))) NE -1) THEN BEGIN
sat='PROBA1'
i_srem_unit=1
ENDIF

IF (((I = STRPOS(PACC_file, 'IREM_'))) NE -1) THEN BEGIN
sat='INTEGRAL'
i_srem_unit=2
ENDIF

IF (((I = STRPOS(PACC_file, 'SREMROSETTA_'))) NE -1) THEN BEGIN
sat='ROSETTA'
i_srem_unit=3
ENDIF

IF (((I = STRPOS(PACC_file, 'SREMHERSCHEL_'))) NE -1) THEN BEGIN
sat='HERSCHEL'
i_srem_unit=4
ENDIF

IF (((I = STRPOS(PACC_file, 'SREMPLANCK_'))) NE -1) THEN BEGIN
sat='PLANCK'
i_srem_unit=5

ENDIF

PRINT, "SATELLITE=",sat
PRINT, "SREM unit=",i_srem_unit

;--------------------------------------------------------------------;
;---------Extract COUNTS from l1_cdf_file_name (PACC file)-----------;
;--------------------------------------------------------------------;
;SET_PADVALUE=0.0
cnts=extract_counts(PACC_file)
help,cnts

;--------------------------------------------------------------------------------------------;
;-----Calculate the FLUXES with the SVD method - this is the CORE of the whole process-------;
;--------------------------------------------------------------------------------------------;
;i_srem_unit=2
fluxes_info=COUNTS_TO_FLUXES(cnts, i_srem_unit)
help,fluxes_info

;------------------------------------------------------------------;
;--------Save the extracted variables from FSVD to SAV files-------;
;------------------------------------------------------------------;
SAVE, fluxes_info.FPDO, fluxes_info.FEDO, fluxes_info.FPDO_ENERGY, fluxes_info.FEDO_ENERGY, fluxes_info.FPDO_ENERGY_BOU, fluxes_info.FEDO_ENERGY_BOU, fluxes_info.PARTICLE_INDEX, fluxes_info.GFP, fluxes_info.GFE, /VARIABLES, FILENAME='/home/sremdc/esa/output/FSVD.sav'
;Restore the SAV file
RESTORE, '/home/sremdc/esa/output/FSVD.sav'


stop



end

