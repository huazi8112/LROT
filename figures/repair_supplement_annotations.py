"""Correct only declared S15 text bands; preserve every data pixel."""
def repair_s15(image):
 from PIL import Image,ImageDraw,ImageFont
 from matplotlib import font_manager
 image=image.copy();font_path=font_manager.findfont('DejaVu Sans');font=ImageFont.truetype(font_path,46);axis_font=ImageFont.truetype(font_path,36)
 boxes=[(185,0,1810,120),(0,100,90,1140)]
 draw=ImageDraw.Draw(image)
 for box in boxes:draw.rectangle((box[0],box[1],box[2]-1,box[3]-1),fill='white')
 draw.text((997,50),'(A) Normalized cost shares at final plan',font=font,fill='black',anchor='mm')
 tile=Image.new('RGB',(850,65),'white');td=ImageDraw.Draw(tile);td.text((425,32),'Normalized weighted-cost share',font=axis_font,fill='black',anchor='mm');tile=tile.rotate(90,expand=True);image.paste(tile,(5,182))
 return image,boxes
