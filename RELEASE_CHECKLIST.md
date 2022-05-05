#Release Checklist

Please be sure to do the following when making a release. This is a template, mostly  it concerns
updates to Readthedocs (RTD), codecov, and Zenodo--all things that have not implemented here yet.

1. (when we have documentation on RTD) Update release page of readthedocs documentation `docs/source/08_Release_History.rst`.
2. (when we have documentation on RTD) Update citation page of readthedocs `docs/source/07_Citations.rst`. (in 3 places--find em' all!)
3. Update the front-facing README in the Version history and citations sections. (in 4 places including version history -- find em' all!)
5. Update release notes `RELEASENOTES.md`.
6. Update `docs/source/index.rst` if badges changed.
7. Restore fail-on-warning on .readthedocs.yaml if it was turned off.
8. Make a release on github.
9. Be sure the `stable` build of readthedocs points to the new release.
10. Be sure to create a version on readthedocs of the new release. 
11. If the version is a patch, deactivate the docs version for the previous patch of the same minor version. (Only one docs version for each minor version should be active at a time.)
12. (when we have codecov) Be sure codecov website is switched to default to master branch.
13. (When we have Zenodo) Update Zenodo.
