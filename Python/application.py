import PAP_Optimal
import PAP_Simple
import numpy as np
import pandas as pd
import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import io # For converting figure to image
import base64

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])
app.title = "The PAP App"
server = app.server

# App layout
app.layout = html.Div([
    html.Div([
        "Source code and help: ", 
        dcc.Link('Github', href='https://github.com/maxkasy/The_PAP_App', target='_blank'),
    ], style={'text-align': 'right'}),
    html.Div([
        "© 2024 ", 
        dcc.Link('Maximilian Kasy', href='https://maxkasy.github.io/home/', target='_blank'),
    ], style={'text-align': 'right'}),
    html.H1('The PAP App', className='text-center mb-4', style={'margin': '50px 0'}),
    dbc.Row([
        dbc.Col(dcc.Dropdown(
            id='data-type',
            options=[
                {'label': 'Binary Data', 'value': 'binary'},
                {'label': 'Normal Data', 'value': 'normal'}
            ],
            placeholder="Select a data type"
        )),
        dbc.Col(dcc.Checklist(id='simple-cutoff', options=[{'label': ' Simple cut-off rule', 'value': 'SCR'}], value=[])),
        dbc.Col(html.Label('Size of test')),
        dbc.Col(dcc.Input(id='input-size', type='number', placeholder="Enter size value", value=.05, step=.05))
    ], className='mb-3'),
    html.Hr(style={'margin': '40px 0'}),
    # Input fields for binary data
    html.Div([
        dbc.Row([
            dbc.Col(html.Label('Prob of observing')),
            dbc.Col(dcc.Input(id='input-etaJ-binary', type='text', placeholder="Enter etaJ values (comma separated)", value=".9,.7,.5")),
            dbc.Col(html.Label('Null hypothesis')),
            dbc.Col(dcc.Input(id='input-minp', type='number', placeholder="Enter minp value", value=.05, step=.05))
        ], className='mb-3'),
        dbc.Row(html.Label('Parameters of Beta prior')),
        dbc.Row([
            dbc.Col(html.Label('alpha')),
            dbc.Col(dcc.Input(id='input-alpha', type='number', placeholder="Enter alpha value", value=1, step=0.1)),
            dbc.Col(html.Label('beta')),
            dbc.Col(dcc.Input(id='input-beta', type='number', placeholder="Enter beta value", value=1, step=0.1))
        ], className='mb-3'),
    ], id='input-fields-binary', style={'display': 'none'}),
    # Input fields for normal data
    html.Div([
        dbc.Row([
            dbc.Col(html.Label('Prob of observing')),
            dbc.Col(dcc.Input(id='input-etaJ-normal', type='text', placeholder="Enter etaJ values (comma separated)", value=".9,.5")),
            dbc.Col(),
            dbc.Col()
        ], className='mb-3'),
        dbc.Row([
            dbc.Col(html.Label('Null mean vector')),
            dbc.Col(dcc.Input(id='input-mu0', type='text', placeholder="Enter mu0 values (comma separated)", value="0,0")),
            dbc.Col(html.Label('Null variance matrix')),
            dbc.Col(dcc.Input(id='input-Sigma0', type='text', placeholder="Enter Sigma0 values (matrix)", value="[[1, 0], [0, 1]]"))
        ], className='mb-3'),
        dbc.Row([
            dbc.Col(html.Label('Prior mean vector')),
            dbc.Col(dcc.Input(id='input-mu', type='text', placeholder="Enter mu values (comma separated)", value="1,1")),
            dbc.Col(html.Label('Prior variance matrix')),
            dbc.Col(dcc.Input(id='input-Sigma', type='text', placeholder="Enter Sigma values (matrix)", value="[[2, 1], [1, 2]]"))
        ], className='mb-3'),
    ], id='input-fields-normal', style={'display': 'none'}),
    dbc.Row([dbc.Col(html.Button('Run', id='run-button')),
        dbc.Col(html.Button("Download PAP", id="btn_csv"))]),
    dcc.Download(id="download-data"),
    html.Hr(style={'margin': '40px 0'}),
    # Output fields
    html.Div([
        html.Div(id='power-binary', style={'text-align': 'center', 'font-size': '20px'}),
        html.Br(),
        html.Div(id='output-binary')
    ], id='output-fields-binary', style={'display': 'none'}),
    html.Div([
        html.Div(id='power-normal', style={'text-align': 'center', 'font-size': '20px'}),
        html.Br(),
        html.Div(id='output-normal'),
        html.Img(id='output-fig-normal', 
         style={'display': 'block', 'marginLeft': 'auto', 'marginRight': 'auto'}),
    ], id='output-fields-normal', style={'display': 'none'}),
    dcc.Store(id='store-t-binary'),
    dcc.Store(id='store-t-normal')
], style={'margin': '10%'})

# Switching input fields based on data type
@app.callback(
    [Output('input-fields-binary', 'style'),
     Output('input-fields-normal', 'style')],
    [Input('data-type', 'value')]
)
def toggle_input_fields(data_type):
    if data_type == 'binary':
        return {'display': 'block'}, {'display': 'none'}
    elif data_type == 'normal':
        return {'display': 'none'}, {'display': 'block'}
    else:
        return {'display': 'none'}, {'display': 'none'}

# Switching output fields based on data type
@app.callback(
    [Output('output-fields-binary', 'style'),
     Output('output-fields-normal', 'style')],
    [Input('data-type', 'value')]
)
def toggle_output(data_type):
    if data_type == 'binary':
        return {'display': 'block'}, {'display': 'none'}
    elif data_type == 'normal':
        return {'display': 'none'}, {'display': 'block'}
    else:
        return {'display': 'none'}, {'display': 'none'}

# Calculating optimal PAP for binary data
@app.callback(
    [Output('power-binary', 'children'),
     Output('output-binary', 'children'),
     Output('store-t-binary', 'data')],
    [Input('run-button', 'n_clicks')],
    [State('data-type', 'value'),
     State('input-etaJ-binary', 'value'),
     State('input-minp', 'value'),
     State('input-alpha', 'value'),
     State('input-beta', 'value'),
     State('input-size', 'value'),
     State('simple-cutoff', 'value')]
)
def update_output_binary(n_clicks, data_type, input_etaJ_binary, input_minp, input_alpha, input_beta, input_size, simple_cutoff):
    if n_clicks is None or data_type != 'binary':
        raise dash.exceptions.PreventUpdate
    else:
        input_binary = {
            'etaJ': np.array([float(i) for i in input_etaJ_binary.split(',')]),
            'minp': float(input_minp),
            'alpha': float(input_alpha),
            'beta': float(input_beta)
        }
        if simple_cutoff == []:
            test_args = PAP_Optimal.pap_binary_data(input_binary)
            opt_test_binary = PAP_Optimal.optimal_test(test_args, input_size)

            power_message = f"Expected power: {opt_test_binary['power'].round(2)}"

            t = opt_test_binary["t"].round(2)
            t_filtered = t.query('t > 0').sort_values('t')                 
        else:
            opt_simple_binary = PAP_Simple.pap_binary_data_simple(input_binary, input_size)
            power_message = f"Expected power: {opt_simple_binary['power'].round(2)}"
            t = opt_simple_binary["test"]  
            t_filtered = t

        table = dash_table.DataTable(
                data=t_filtered.to_dict('records'),
                columns=[{'name': i, 'id': i} for i in t_filtered.columns],
                style_cell={'minWidth': '100px', 'maxWidth': '200px'},
                style_table={'maxHeight': '300px', 'overflowY': 'auto'},
            )  
         
        return power_message, table, t.to_csv(index=False)

# Calculating optimal PAP for normal data
@app.callback(
    [Output('power-normal', 'children'),
     Output('output-normal', 'children'),
     Output('store-t-normal', 'data'),
     Output('output-fig-normal', 'src')],
    [Input('run-button', 'n_clicks')],
    [State('data-type', 'value'),
     State('input-etaJ-normal', 'value'),
     State('input-mu0', 'value'),
     State('input-Sigma0', 'value'),
     State('input-mu', 'value'),
     State('input-Sigma', 'value'),
     State('input-size', 'value'),
     State('simple-cutoff', 'value')]
)
def update_output_normal(n_clicks, data_type, input_etaJ_normal, input_mu0, input_Sigma0, input_mu, input_Sigma, input_size, simple_cutoff):
    img = dash.no_update  # Initialize img to dash.no_update

    if n_clicks is None or data_type != 'normal':
        raise dash.exceptions.PreventUpdate
    else:
        input_normal = {
            'etaJ': np.array([float(i) for i in input_etaJ_normal.split(',')], ndmin=1),
            'mu0': np.array([float(i) for i in input_mu0.split(',')], ndmin=1),
            'Sigma0': np.array(eval(input_Sigma0), ndmin=2),
            'mu': np.array([float(i) for i in input_mu.split(',')], ndmin=1),
            'Sigma': np.array(eval(input_Sigma), ndmin=2)
        }
        if simple_cutoff == []:
            test_args = PAP_Optimal.pap_normal_data(input_normal)
            opt_test_normal = PAP_Optimal.optimal_test(test_args, input_size)

            power_message = f"Expected power: {opt_test_normal['power'].round(2)}"

            t = opt_test_normal["t"].round(2)
            t_filtered = t.query('t > 0').sort_values('t')     

            if len(input_normal["etaJ"]) == 2:
                fig = PAP_Optimal.plot_normal_pap(t)
                # Convert the matplotlib figure to a base64 image
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                image_base64 = base64.b64encode(buf.read()).decode()
                img = 'data:image/png;base64,{}'.format(image_base64)
        else:
            opt_simple_normal = PAP_Simple.pap_normal_data_simple(input_normal, input_size)
            power_message = f"Expected power: {opt_simple_normal['power'].round(2)}"
            t = opt_simple_normal["test"]     
            t_filtered = t
        
        table = dash_table.DataTable(
                data=t_filtered.to_dict('records'),
                columns=[{'name': i, 'id': i} for i in t_filtered.columns],
                style_cell={'minWidth': '100px', 'maxWidth': '200px'},
                style_table={'maxHeight': '300px', 'overflowY': 'auto'}
            )  
        
        return power_message, table, t.to_csv(index=False), img


# Download PAP
@app.callback(
    Output("download-data", "data"),
    [Input("btn_csv", "n_clicks")],
    [State('data-type', 'value'),
     State('store-t-binary', 'data'),
     State('store-t-normal', 'data')],
    prevent_initial_call=True,
)
def download_pap(n_clicks, data_type, t_binary, t_normal):
    if n_clicks is None:
        raise dash.exceptions.PreventUpdate

    if data_type == 'binary':
        return dict(content=t_binary, filename="pap_binary.csv", type='text/csv')
    elif data_type == 'normal':
        return dict(content=t_normal, filename="pap_normal.csv", type='text/csv')



if __name__ == '__main__':
    app.run_server(debug=True)