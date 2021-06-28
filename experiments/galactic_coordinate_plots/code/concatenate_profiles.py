from PIL import Image, ImageOps, ImageDraw, ImageFont
path = "output_heliocentric_profile_5_kpc/clean_streams_reduced_heliocentric_distance_"
# the paths
image_list = ["0_5", "5_10", "10_15", "15_20", "20_25", "25_30", "30_35", "45_300"]
imgs    = [ Image.open(path+i+"_kpc.png") for i in image_list ]
# define the size
img_w,img_h=imgs[0].size
rows = 4
cols = 2
# crop
left, img_h = imgs[0].crop((150,0,.905*img_w, img_h)).size
right,img_h = imgs[1].crop((290,0,img_w, img_h)).size
# define shape of the final image
grid_w = left+right 
grid_h = img_h*rows
# create new image
grid = Image.new("RGB", size=(grid_w, grid_h), color=(255,255,255))
grid_w, grid_h = grid.size
# paste the sub images
for i, img in enumerate(imgs):
    if i%cols==0:
        grid.paste(img.crop((150,0,.905*img_w, img_h)), box=(0, i//cols*img_h))
    else:
        grid.paste(img.crop((290,0,img_w, img_h)), box=(grid_w//2 , i//cols*img_h))
grid.save("../output_heliocentric_profile/heliocentric_slices.png")