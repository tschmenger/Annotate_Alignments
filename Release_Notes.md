# V11 (04. October 2024)
- Fixed feature annotation text incorrectly overlapping when there are more than 9 different features
- added % conservation mouseover text to the black sequence conservation bars

# V8 (09. August 2023)
Applied several changes to the script.
- added progress notifications
- the script returns dictionaries (for positional and structural information) if none were provided by the user. This enables the user to re-use them and easily try different settings without having to fetch the information from the web again
- removed the requirement of the 'svgwrite' module by using raw svg code instead
- removed the requirement of the "Biopython" module by parsing the alignmentfile manually
