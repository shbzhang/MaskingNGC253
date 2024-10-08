This is a simple script to get mask of emission for NGC253.

We have datacube for NGC253 in H13CN(1-0), HCNH+(2-1), (3-2), (4-3), and (5-4), but some velocity components are contaminated by other emission lines within the band along the line of sight.

To solve this, we use HCNH+(2-1) as an initial mask since it has no other strong lines in its band. The mask is then expanded and reprojected to other observed lines.

HCNH+(4-3) needs further attention. Strong contaminations are found at the GMC-3, and 5. So their masks are further removed from the mask generated in last step.
