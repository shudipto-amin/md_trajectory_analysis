# README

Run script `analyse.py` with no arguments to calculate RDF between oxygen-oxygen in the `example.pdb`.
This will generate a plot of RDF, $$g_{ab}$$, and coordination number, $$N_{ab}$$ as function of $$r$$ in `rdf_and_coordination.png`.

Often by 'coordination number' people actually mean the value of $$N_{ab}$$ at $$r$$ where first peak of $$g_{ab}$$ occurs. 
This you can either infer from the plot, or write code find the first maxima $$r_{max}$$ and then get $$N_{ab}(r_{max})$$. 
Similarly, for second coordination shell, you can use the second peak and so on. 
