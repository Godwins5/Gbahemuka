 classdef (Abstract) NumericalMethod
    % Abstract parent class def for mumerical methods
    properties (Access = protected)
        f  % fuction Handle
        a  % lower limit
        b  % upper limit
        n  % number of steps/interval
    end

    methods
        
       function obj = Numericalmethods(f,a,b,n)
            obj.f = f;
            obj.a = a;
            obj.b = b;
            obj.n = n;
        end
    end

    methods (Abstract)
        result = compute(obj);
    end
 end

 %% Gui for differential solver 
 classdef DifferentialSolver < NumericalMethod
    % solves dy/dx = f(x,y) using eulers method

    properties
        y0 % initial condition
    end

    methods
        function obj = DifferentialSolver(f,a,b,n,y0)
            obj@NumericalMethod(f,a,b,n);
            obj.y0 = y0;
        end

        function y = compute(obj)
            h = (obj.b - obj.a);
            x = obj.a:h:obj.b;
            y = zeros(1,obj.n + 1);
            y(1) = obj.y0;

            for i = 1 : obj.n
                y(i +1) = y(i) + h*obj.f(x(i),y(i));
            end

            plot(x,y,"-o");
            title('Euler method solution');
            xlabel("x");
            ylabel("y");
            grid on;
        end
    end
end

%% gui codes


    

function de_solver_gui()
    % MATLAB GUI for solving ODEs using Euler's method
    % User inputs: f(t,y) expression, y0, t0, tf, h
    % Outputs: Plot of y vs t

    % Create main figure
    fig = figure('Name', 'Euler Method ODE Solver', 'NumberTitle', 'off', ...
                 'Position', [100, 100, 500, 400], 'MenuBar', 'none');

    % Axes for plotting
    ax = axes('Parent', fig, 'Position', [0.1, 0.2, 0.8, 0.6]);

    % UI Controls
    % f(t,y) input
    uicontrol('Style', 'text', 'Position', [10, 350, 80, 20], 'String', 'f(t,y):', ...
              'Parent', fig, 'FontSize', 10);
    edit_f = uicontrol('Style', 'edit', 'Position', [90, 350, 300, 20], ...
                       'String', 'sin(t) - 2*y', 'Parent', fig);  % Example: y' = sin(t) - 2y

    % y0 input
    uicontrol('Style', 'text', 'Position', [10, 320, 80, 20], 'String', 'y0:', ...
              'Parent', fig, 'FontSize', 10);
    edit_y0 = uicontrol('Style', 'edit', 'Position', [90, 320, 50, 20], ...
                        'String', '1', 'Parent', fig);

    % t0 input
    uicontrol('Style', 'text', 'Position', [160, 320, 80, 20], 'String', 't0:', ...
              'Parent', fig, 'FontSize', 10);
    edit_t0 = uicontrol('Style', 'edit', 'Position', [200, 320, 50, 20], ...
                        'String', '0', 'Parent', fig);

    % tf input
    uicontrol('Style', 'text', 'Position', [270, 320, 80, 20], 'String', 'tf:', ...
              'Parent', fig, 'FontSize', 10);
    edit_tf = uicontrol('Style', 'edit', 'Position', [300, 320, 50, 20], ...
                        'String', '5', 'Parent', fig);

    % h input
    uicontrol('Style', 'text', 'Position', [370, 320, 80, 20], 'String', 'h:', ...
              'Parent', fig, 'FontSize', 10);
    edit_h = uicontrol('Style', 'edit', 'Position', [400, 320, 50, 20], ...
                       'String', '0.1', 'Parent', fig);

    % Solve button
    btn_solve = uicontrol('Style', 'pushbutton', 'Position', [200, 280, 80, 30], ...
                          'String', 'Solve & Plot', 'Parent', fig, 'Callback', @solve_callback);

    % Callback function
    function solve_callback(~, ~)
        % Get input values
        f_expr = get(edit_f, 'String');
        y0 = str2double(get(edit_y0, 'String'));
        t0 = str2double(get(edit_t0, 'String'));
        tf = str2double(get(edit_tf, 'String'));
        h = str2double(get(edit_h, 'String'));

        % Validate inputs
        if isnan(y0) || isnan(t0) || isnan(tf) || isnan(h) || h <= 0
            errordlg('Invalid input values. Ensure y0, t0, tf, h are numbers and h > 0.', 'Input Error');
            return;
        end
        if t0 >= tf
            errordlg('t0 must be less than tf.', 'Input Error');
            return;
        end

        % Create anonymous function for f(t,y)
        try
            f = str2func(['@(t,y) ' f_expr]);
        catch
            errordlg('Invalid f(t,y) expression. Use MATLAB syntax, e.g., sin(t) - y.', 'Syntax Error');
            return;
        end

        % Euler's method implementation
        n_steps = floor((tf - t0) / h) + 1;
        t = linspace(t0, tf, n_steps);
        y = zeros(size(t));
        y(1) = y0;
        for i = 1:n_steps-1
            y(i+1) = y(i) + h * f(t(i), y(i));
        end

        % Plot
        cla(ax);  % Clear axes
        plot(ax, t, y, 'b-', 'LineWidth', 2);
        xlabel(ax, 't');
        ylabel(ax, 'y');
        title(ax, 'Solution using Euler''s Method');
        grid(ax, 'on');
    end
end

%% intergral solver
classdef IntergralSolver < NumericalMethod
    % computes intergral using Trapezoidal Rule

    methods
        function result = compute(obj)
            h = (obj.b - obj.a)/obj.n;
            x = obj.a:h:obj.b;
            y = obj.f(x);
            result = (h/2)*(y(1) + 2*sum(y(2:end - 1)) + y(end));
            
            fprintf("Approximate Intergral value = %.5f\n",result);
        end
    end
end

%% gui for intergral solver
function TrapezoidalIntegratorApp()
    % Create the main figure window
    fig = uifigure('Name', 'Trapezoidal Rule Integrator', ...
                   'Position', [100, 100, 950, 700], ...
                   'Color', [0.94, 0.94, 0.96]);

    % Header label
    header = uilabel(fig, 'Text', 'Definite Integral using Trapezoidal Rule', ...
                     'FontSize', 18, 'FontWeight', 'bold', ...
                     'Position', [50, 650, 850, 40], ...
                     'HorizontalAlignment', 'center');

    % Input panel
    inputPanel = uipanel(fig, 'Title', 'Inputs', ...
                         'Position', [50, 400, 850, 200], ...
                         'BackgroundColor', [0.94, 0.94, 0.96]);

    % Function input
    uilabel(inputPanel, 'Text', 'f(x) =', 'Position', [20, 150, 50, 20], ...
            'FontSize', 12);
    funcEdit = uieditfield(inputPanel, 'text', 'Value', 'x.^2 + sin(x)', ...
                           'Position', [80, 150, 300, 20], ...
                           'FontName', 'Courier New');

    % Lower limit a
    uilabel(inputPanel, 'Text', 'a =', 'Position', [20, 120, 30, 20], ...
            'FontSize', 12);
    aEdit = uieditfield(inputPanel, 'numeric', 'Value', 0, ...
                        'Position', [60, 120, 80, 20], ...
                        'FontSize', 12);

    % Upper limit b
    uilabel(inputPanel, 'Text', 'b =', 'Position', [160, 120, 30, 20], ...
            'FontSize', 12);
    bEdit = uieditfield(inputPanel, 'numeric', 'Value', 2, ...
                        'Position', [200, 120, 80, 20], ...
                        'FontSize', 12);

    % Number of trapezoids n
    uilabel(inputPanel, 'Text', 'n =', 'Position', [20, 90, 30, 20], ...
            'FontSize', 12);
    nEdit = uieditfield(inputPanel, 'numeric', 'Value', 100, ...
                        'Position', [60, 90, 80, 20], ...
                        'FontSize', 12);

    % Compute button
    computeBtn = uibutton(inputPanel, 'push', 'Text', 'Compute Integral', ...
                          'Position', [350, 50, 150, 30], ...
                          'FontSize', 12, 'FontWeight', 'bold', ...
                          'BackgroundColor', [0.29, 0.56, 0.89], ...
                          'ButtonPushedFcn', @(btn,event) computeIntegral());

    % Results panel
    resultsPanel = uipanel(fig, 'Title', 'Results', ...
                           'Position', [50, 150, 850, 200], ...
                           'BackgroundColor', [0.94, 0.94, 0.96]);

    resultsText = uitextarea(resultsPanel, 'Value', {'Enter values and click Compute.'}, ...
                             'Position', [20, 20, 810, 160], ...
                             'FontSize', 11, 'Editable', 'off');

    % Plot axes
    ax = uiaxes(fig, 'Position', [50, 20, 850, 120]);
    title(ax, 'Function Plot');
    xlabel(ax, 'x');
    ylabel(ax, 'f(x)');
    grid(ax, 'on');

    % Trapezoidal rule computation function
    function computeIntegral()
        try
            % Parse inputs
            funcStr = funcEdit.Value;
            a = aEdit.Value;
            b = bEdit.Value;
            n = round(nEdit.Value); % Ensure integer
            if n < 1
                n = 1;
            end

            % Create symbolic function for exact integral (optional, for comparison)
            syms x;
            f_sym = str2sym(funcStr);
            exactIntegral = double(vpa(int(f_sym, a, b)));

            % Numerical trapezoidal rule
            h = (b - a) / n;
            x_vals = a:h:b;
            y_vals = double(subs(f_sym, x, x_vals));

            % Trapezoidal approximation
            trapIntegral = h * (0.5 * y_vals(1) + sum(y_vals(2:end-1)) + 0.5 * y_vals(end));

            % Error
            error = abs(trapIntegral - exactIntegral);

            % Update results text
            resultsStr = {sprintf('Approximate Integral: %.6f', trapIntegral), ...
                          sprintf('Exact Integral: %.6f', exactIntegral), ...
                          sprintf('Absolute Error: %.2e', error), ...
                          sprintf('h = %.4f, n = %d', h, n)};
            resultsText.Value = resultsStr;

            % Plot the function and trapezoids
            cla(ax);
            hold(ax, 'on');
            plot(ax, x_vals, y_vals, 'b-', 'LineWidth', 2); % Function curve

            % Plot trapezoids
            for i = 1:n
                x_trap = [x_vals(i), x_vals(i), x_vals(i+1), x_vals(i+1)];
                y_trap = [y_vals(i), y_vals(i+1), y_vals(i+1), y_vals(i)];
                fill(ax, x_trap, y_trap, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
            end

            hold(ax, 'off');
            xlim(ax, [a, b]);
            ylim(ax, 'auto');

        catch ME
            resultsText.Value = {'Error: ' + ME.message};
        end
    end
end


