
export SPHINX_APIDOC_OPTIONS='members,private-members,undoc-members,show-inheritance'

echo `rm -f ./sphinx/*.html`
echo `rm -f ./sphinx/searchindex.js ./sphinx/objects.inv`
echo `rm -f ./sphinx/src/*.rst`
echo `rm -fR ./sphinx/_modules/*`
echo `rm -fR ./sphinx/_sources/*`
echo `rm -fR ./sphinx/_static/*`

echo `sphinx-apidoc -o ./sphinx/src ../python/mangadap --full --separate -H mangadap -A "SDSS-IV/MaNGA Pipeline Group" -V 2.2.2 -R 2.2.2`

echo `sphinx-build -aEb html -w ./sphinx/warnings.out -T ./sphinx/src ./sphinx`

