
# Distribution of number of atoms around Zn within cutoff
cutoff = 2.1  # choose based on RDF
atom = "ZN"
seltxt = f"around {cutoff} name {atom}"
counts, bin_edges = uni.get_contact(seltxt)

fig, ax = pp.subplots()
ax.set_xlabel("Coordination number")
ax.set_ylabel("Frequency")
ax.hist(counts, bin_edges)
if save:
    fig.savefig(f"contact_around_{atom}_{cutoff}.png")

fig,ax = pp.subplots()
ax.set_xlabel("Frame no.")
ax.set_ylabel("Coordination number")
ax.plot(counts)
if save:
    fig.savefig(f"cn_vs_frame_{atom}_{cutoff}.png")

