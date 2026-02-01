from PIL import Image

folder = "circle_data_u=1"
file = "/friction_az"

ls = []
for i in range(50):
    name = (
        f"/Users/hdknkms/Desktop/地球惑星科学課題演習DC/last_presentation/"
        + folder
        + file
        + f"_{i*100}.png"
    )
    im = Image.open(name)
    ls.append(im)

ls[0].save(
    "/Users/hdknkms/Desktop/地球惑星科学課題演習DC/last_presentation/"
    + folder
    + file
    + ".gif",
    save_all=True,
    append_images=ls[1:],
    optimize=False,
    duration=300,
    loop=0,
)
