def save_figure(fig, fn, width=800, height=600, scale=2):
    fig.write_image(f"{fn}.png")
    fig.write_image(f"{fn}.pdf", format="pdf")