This is a simple script to get mask of emission for NGC253.

We have datacubes for NGC253 observed in several molecular lines: H13CN(1-0), HCNH+(2-1), HCNH+(3-2), HCNH+(4-3), and HCNH+(5-4). However, some velocity components are contaminated by other emission lines within the band along the line of sight, making it challenging to isolate the emission from each transition.

To address this, we use the HCNH+(2-1) line, which is free from strong contaminating lines, as an initial mask. This mask is then expanded and reprojected onto the other observed lines.

For the HCNH+(4-3) transition, further adjustments are required due to significant contamination in regions such as GMC-3 and GMC-5. The emission in these regions is manually excluded from the mask generated in the previous step.
