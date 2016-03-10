	

module funcapp

	export runall


	using PyPlot
	using ApproxFun
	using FastGaussQuadrature
	# using ApproXD

	function chebBasis(grid, n, lb, ub)
		Phi = zeros(length(grid), n)
		for i = 1:length(grid)
			# x = 2*(grid[i]-lb)/(ub-lb)-1
			for j = 1:n
				Phi[i,j] = cos(acos(grid[i])*(j-1))
				# Phi[i,j] = cos((n-i+0.5)*(j-1)*pi/n)
			end
		end
		return Phi
	end

	

	# use chebyshev to interpolate this:
	function q1(n)

		f(x) = x + 2*x.^2 - exp(-x)	
		lb 		= -3;
		ub 		= 3;
		grid 	= linspace(lb, ub, n)
		ngrid 	= length(grid)
		trueFun = f(grid)

		# nods 			= gausschebyshev(n)
		nodes = zeros(n,1)
		for k = 1:n
			nodes[k] 			= cos(((2*k-1)/2n)*pi)
		end
	
		rescaleNodes 	= (1/2)*(ub+lb)+(1/2)*(ub-lb)*nodes		
		Phi = chebBasis(nodes, n, lb, ub)

		Y = f(rescaleNodes)
		chebCoeffs = inv(Phi)*Y
		rescaleGrid =  (1/2)*(ub+lb)+(1/2)*(ub-lb)*linspace(-1,1,50)	
		
		global interpPhi = chebBasis(linspace(-1,1,50), n, lb, ub)
		approxY = interpPhi*chebCoeffs
		newY = f(rescaleGrid)
		errorY = newY - approxY
		figure()
		subplot(211)
		plot([newY approxY])
		subplot(212)
		plot(errorY)
		title("Q1 Error")

	end

	function q2(n)
		f(x) = x + 2*x.^2 - exp(-x)	
		lb 		= -3;
		ub 		= 3;
		d = Space([-3,3]) ;
		pts = points(d, n);
		y = f(pts);
		chebApprox = Fun(ApproxFun.transform(d, y), d)
		grid 	= linspace(lb, ub, n)
		evalap=chebApprox(collect(grid))
		error=abs(evalap-f(grid))
		figure()
		plot(grid,error)
		title("Q2 Error")

	end


	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()
	
		figure()
		subplot(331)
		plot(interpPhi[:,1])
		subplot(332)
		plot(interpPhi[:,2])
		subplot(333)
		plot(interpPhi[:,3])
		subplot(334)
		plot(interpPhi[:,4])
		subplot(335)
		plot(interpPhi[:,5])
		subplot(336)
		plot(interpPhi[:,6])
		subplot(337)
		plot(interpPhi[:,7])
		subplot(338)
		plot(interpPhi[:,8])
		subplot(339)
		plot(interpPhi[:,9])


	end

	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	type ChebyType
		f::Function # fuction to approximate 
		nodes::Union{Vector,LinSpace} # evaluation points
		basis::Matrix # basis evaluated at nodes
		coefs::Vector # estimated coefficients

		deg::Int 	# degree of chebypolynomial
		lb::Float64 # bounds
		ub::Float64

		# constructor
		function ChebyType(_nodes::Union{Vector,LinSpace},_deg,_lb,_ub,_f::Function)
			n = length(_nodes)
			y = _f(_nodes)
			_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
			_coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
			# create a ChebyType with those values
			new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
		end
	end
	
	# function to predict points using info stored in ChebyType
	function predict(Ch::ChebyType,x_new)

		true_new = Ch.f(x_new)
		basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds = basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end

	Runge(x) = 1./(1+25*x.^2)
	# function q4a(deg=(5,9,15),lb=-1.0,ub=1.0)
	function q4a()
		pts = 60
		Y = zeros(pts, 3)
		newY = Runge(linspace(-5,5,pts))

		approxChebY = zeros(pts, 3)
		approxUniY = zeros(pts, 3)
		gDeg = [5, 9, 15]
		for i  = 1:3
			
			deg = gDeg[i]
			uniformNodes 		= linspace(-1,1, deg)
			chebNodes = zeros(deg, 1)
			for k = 1:deg
				chebNodes[k] 			= cos(((2*k-1)*pi/(2*deg)))
			end
			lb 		= -5;
			ub 		= 5;
			# chebNodes    		= gausschebyshev(5)
			rescaleChebNodes 	=  (1/2)*(ub+lb)+(1/2)*(ub-lb)*chebNodes
			rescaleUniNodes 	=  (1/2)*(ub+lb)+(1/2)*(ub-lb)*uniformNodes

		
			B 		= chebBasis(chebNodes, deg, lb, ub)
			uniB 	= chebBasis(uniformNodes, deg, lb, ub)
			chebY = Runge(rescaleChebNodes)
			uniY = Runge(rescaleUniNodes)
			chebCoeffs 		= inv(B)*chebY
			uniCoeffs 		= inv(uniB)*uniY
			interpB 		= chebBasis(linspace(-1,1,pts), deg, lb, ub)
			approxChebY[:,i] 	= interpB*chebCoeffs
			approxUniY[:,i] 		= interpB*uniCoeffs
			
		end
		figure()
		subplot(211)
		plot(linspace(-5,5,pts), [newY approxChebY])
		subplot(212)
		plot(linspace(-5,5,pts),[newY approxUniY])

	end

	function q4b()
		
		# bs = BSpline(13,3,-5,5)
		# B = full(getBasis(collect(linspace(-5,5.0,500)),bs))
		println("q4 not completed")
		
	end

	function q5()
		println("q5 not completed")

		
	end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		q1(15)
		q2(15)
		q3()
		q4a()
		q4b()
		q5()
	end


end

