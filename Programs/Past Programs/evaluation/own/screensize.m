function data = screensize(dx,dy,px,py)

data.pixel_dxmum = dx/px*1e3;
data.pixel_dymum = dy/py*1e3;
data.diagonal_mm = sqrt(dx^2+dy^2);
data.diagonal_in = data.diagonal_mm/25.4;

end