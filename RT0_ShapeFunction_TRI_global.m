clc
clear all
close all

% coordinates of the source element
coord = [    3.2217    9.8066
    4.6208    8.1898
    4.9911   10.3278];

figure(1)

title("Source triangle")
patch('Vertices', coord, 'Faces', 1:3, 'FaceVertexCData', 0, 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
xlabel('x')
ylabel('y')
pbaspect([1, 1, 1])
for i = 1:3
    text(coord(i, 1), coord(i, 2), ['(', num2str(i), ')'], 'color', 'r'); hold on
end

element = [1, 2, 3];


