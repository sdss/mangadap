;+
; NAME:
;   mangadap_version
; PURPOSE:
;   Return the version name for the mangadap product
; CALLING SEQUENCE:
;   vers = mangadap_version()
; OPTIONAL INPUTS:
;   /SIMPLE: Don't return svn revision number if version is trunk
; OUTPUTS:
;   vers       - Version name for the product mangadap
; COMMENTS:
;   Depends on shell script in $MANGADAP_DIR/bin
;-
;------------------------------------------------------------------------------
function mangadap_version,simple=simple
   cmd = "mangadap_version"
   spawn, cmd, stdout, /noshell
   mangadap_version = stdout[0]
   if keyword_set(simple) then mangadap_version=(strsplit(mangadap_version,/extract,' '))[0]
   return, mangadap_version
end
;------------------------------------------------------------------------------
