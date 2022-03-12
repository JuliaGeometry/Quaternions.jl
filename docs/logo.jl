using Luxor, Colors

L = 500
R = 130
r = 25
t = 18
h = 0.33

path_svg = joinpath(@__DIR__, "src", "assets", "logo.svg")
path_ico = joinpath(@__DIR__, "src", "assets", "favicon.ico")

@svg begin
	Drawing(L, L, path_svg)
	setlinecap("round")
	setlinejoin("round")
	p₊ = Point(cos(π/6), -sin(π/6)) * 2R/√3
	p₋ = -p₊/2
	angles = (-2π/3, -π/3, 0, π/3, 4π/3)

	for i in 1:3
		origin(L/2, L/2)
		rotate(4π/3*(i-1))
		sethue(Colors.JULIA_LOGO_COLORS[(:red, :purple, :green)[i]])
		circle(p₋, r, :fill)
		circle(p₊, r, :fill)
		setline(t)
		for j in 1:4
			arc(p₋, R, angles[j]+h, angles[j+1]-h, :stroke)
		end
		setline(t/2)
		poly(ngon(p₋*(1+√3), t, 3, π), :stroke, close=true)
	end
end

run(`convert -density 256x256 -background transparent $(path_svg) -define icon:auto-resize -colors 256 $(path_ico)`)
