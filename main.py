import sys

from sympy import Symbol, lambdify, solve
from sympy.parsing.sympy_parser import parse_expr

from Algorithms import *
from presentation import prepare_contour, summary_to_graph, summary_simplex_to_graph, gamma_to_graph

if __name__ == '__main__':

    LSP = '1813056'
    LSP_a = int(LSP[5])
    LSP_b = int(LSP[6])

    variables = {'x': Symbol('x', real=True),
                 'y': Symbol('y', real=True),
                 'z': Symbol('z', real=True)
                 }
    variables_f = {
        'fx': Symbol('fx', real=True),
        'g1x': Symbol('g1x', real=True),
        'h1x': Symbol('h1x', real=True),
        'h2x': Symbol('h2x', real=True),
        'h3x': Symbol('h3x', real=True)
    }
    form_fx = parse_expr('-(x * y * z)', variables)
    form_g1x = parse_expr('2 * (x * y + x * z + y * z) - 1', variables)
    form_h1x = parse_expr('x * y * -1', variables)
    form_h2x = parse_expr('x * z * -1', variables)
    form_h3x = parse_expr('x * y * z * -1', variables)
    form_bx = parse_expr('g1x**2 ')

    # form_z = solve(form_fx_restriction, variables['z'])[0]
    # form_fx = form_fx.subs(variables['z'], form_z)
    # variables.pop('z')

    # fx = lambdify(variables.values(), form_fx, "numpy")
    # form_dfx = [form_fx.diff(the_symbol) for the_symbol in variables.values()]
    # dfx = [lambdify([v for v in variables.values()],
    #                 the_dfx,
    #                 "numpy")
    #        for the_dfx in form_dfx]
    #
    # def gradient_fx(values):
    #     return np.array([de_fx(*values) for de_fx in dfx])
    #
    # def two_arguments_to_edges(the_x, the_y):
    #     the_z = form_z.subs([(variables['x'], the_x), (variables['y'], the_y)])
    #     the_mult = the_z * the_x * the_y
    #     if the_mult == 0:
    #         return 0, 0, 0
    #     the_length = (the_x * the_z / the_y / 2)**0.5
    #     the_width = (the_y * the_z / the_x / 2)**0.5
    #     the_height = (the_x * the_y / the_z / 2)**0.5
    #     return the_length, the_width, the_height
    #
    # print(f'LSP: {LSP}')
    # print(f'a={LSP_a} ; b={LSP_b}')
    # print('\n\n')
    #
    # near_zero = sys.float_info.min * sys.float_info.epsilon
    # near_zero2 = 1e-16
    # near_zero3 = 1e-6
    # start_length = 0.09
    # gamma_step = 3.8  # 3.8
    # gradient_stop = near_zero2
    # stop_line = 1e-8  # 1e-8
    #
    # x0 = [0, 0]
    # x1 = [1, 1]
    # xm = [LSP_a/10, LSP_b/10]
    # xs = [x0, x1, xm]
    #
    # print("X, f(X), gradf(X)")
    # for x in xs:
    #     print(f"{x}, "
    #           f"{fx(*x)}, "
    #           f"{gradient_fx(x)}")
    #
    # xs_summary = []
    # for x in xs:
    #     xs_summary.append(
    #         [
    #             gradient_descend(fx, gradient_fx, x, stop_line, gamma_step=gamma_step),
    #             the_fastest_descend(fx, gradient_fx, x, stop_line),
    #             deformed_simplex(fx, x, start_length, stop_line, step_limit=189999)
    #         ]
    #     )
    #
    # x_experiments = [(xs[ii], xs_summary[ii]) for ii in range(len(xs))]
    #
    # for start_point, summaries in x_experiments:
    #     for summary in summaries:
    #         summary.translated = two_arguments_to_edges(*summary.solution)
    #         volume = 1
    #         for edge in summary.translated:
    #             volume *= edge
    #         summary.translated_fx = volume
    #     print(f"\n-------------Summary---------\n"
    #           f"Starting point = {start_point}")
    #     print_summary(*summaries)
    #     print(f'-----------------------------\n')
    #
