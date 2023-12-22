import PAP_Optimal
import PAP_Simple
import numpy as np
import pandas as pd
import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])
app.title = "The PAP App"
server = app.server

app.layout = html.Div([
    html.H1('The PAP App', className='text-center mb-4'),
    dcc.Dropdown(
        id='data-type',
        options=[
            {'label': 'Binary Data', 'value': 'binary'},
            {'label': 'Normal Data', 'value': 'normal'}
        ],
        placeholder="Select a data type"
    ),
    html.Hr(),
    html.Div([
        dbc.Row([
            dbc.Col(html.Label('Prob of observing')),
            dbc.Col(dcc.Input(id='input-etaJ-binary', type='text', placeholder="Enter etaJ values (comma separated)", value=".9,.7,.5")),
            dbc.Col(html.Label('Null hypothesis')),
            dbc.Col(dcc.Input(id='input-minp', type='number', placeholder="Enter minp value", value=".05", step="0.05"))
        ], className='mb-3'),
        dbc.Row(html.Label('Parameters of Beta prior')),
        dbc.Row([
            dbc.Col(html.Label('alpha')),
            dbc.Col(dcc.Input(id='input-alpha', type='number', placeholder="Enter alpha value", value="1", step="0.1")),
            dbc.Col(html.Label('beta')),
            dbc.Col(dcc.Input(id='input-beta', type='number', placeholder="Enter beta value", value="1", step="0.1"))
        ], className='mb-3'),
        dbc.Row([dcc.Checklist(id='simple-cutoff-binary', options=[{'label': ' Simple cut-off rule', 'value': 'SCR'}], value=[])],  className='mb-3'),
    ], id='input-fields-binary', style={'display': 'none'}),
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
        dbc.Row([dcc.Checklist(id='simple-cutoff-normal', options=[{'label': ' Simple cut-off rule', 'value': 'SCR'}], value=[])],  className='mb-3'),
    ], id='input-fields-normal', style={'display': 'none'}),
    dbc.Row([dbc.Col(html.Button('Run', id='run-button')),
        dbc.Col(html.Button("Download PAP", id="btn_csv", style={'margin-bottom': '20px'}))]),
    dcc.Download(id="download-data"),
    html.Hr(),
    html.Div([
        html.Div(id='power-binary', style={'text-align': 'center', 'font-size': '20px'}),
        html.Br(),
        html.Div(id='output-binary')
    ], id='output-fields-binary', style={'display': 'none'}),
    html.Div([
        html.Div(id='power-normal', style={'text-align': 'center', 'font-size': '20px'}),
        html.Br(),
        html.Div(id='output-normal')
    ], id='output-fields-normal', style={'display': 'none'}),   
    dcc.Store(id='store-t-binary'),
    dcc.Store(id='store-t-normal')
], style={'margin': '10%'})

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
     State('simple-cutoff-binary', 'value')]
)
def update_output_binary(n_clicks, data_type, input_etaJ_binary, input_minp, input_alpha, input_beta, simple_cutoff_binary):
    if n_clicks is None or data_type != 'binary':
        raise dash.exceptions.PreventUpdate
    else:
        input_binary = {
            'etaJ': np.array([float(i) for i in input_etaJ_binary.split(',')]),
            'minp': float(input_minp),
            'alpha': float(input_alpha),
            'beta': float(input_beta)
        }
        if simple_cutoff_binary == []:
            test_args = PAP_Optimal.pap_binary_data(input_binary)
            opt_test_binary = PAP_Optimal.optimal_test(test_args, .05)

            power_message = f"Expected power: {opt_test_binary['power'].round(2)}"

            t = opt_test_binary["t"].round(2)
            t_filtered = t.query('t > 0').sort_values('t')                 
        else:
            opt_simple_binary = PAP_Simple.pap_binary_data_simple(input_binary, .05)
            power_message = f"Expected power: {opt_simple_binary['power'].round(2)}"
            t = opt_simple_binary["test"]  
            t_filtered = t

        table = dash_table.DataTable(
                data=t_filtered.to_dict('records'),
                columns=[{'name': i, 'id': i} for i in t_filtered.columns],
                style_cell={'minWidth': '100px', 'maxWidth': '200px'},
            )  
         
        return power_message, table, t.to_csv(index=False)





@app.callback(
    [Output('power-normal', 'children'),
     Output('output-normal', 'children'),
     Output('store-t-normal', 'data')],
    [Input('run-button', 'n_clicks')],
    [State('data-type', 'value'),
     State('input-etaJ-normal', 'value'),
     State('input-mu0', 'value'),
     State('input-Sigma0', 'value'),
     State('input-mu', 'value'),
     State('input-Sigma', 'value'),
     State('simple-cutoff-normal', 'value')]
)
def update_output_normal(n_clicks, data_type, input_etaJ_normal, input_mu0, input_Sigma0, input_mu, input_Sigma, simple_cutoff_normal):
    if n_clicks is None or data_type != 'normal':
        raise dash.exceptions.PreventUpdate
    else:
        input_normal = {
            'etaJ': np.array([float(i) for i in input_etaJ_normal.split(',')]),
            'mu0': np.array([float(i) for i in input_mu0.split(',')]),
            'Sigma0': np.array(eval(input_Sigma0)),
            'mu': np.array([float(i) for i in input_mu.split(',')]),
            'Sigma': np.array(eval(input_Sigma))
        }
        if simple_cutoff_normal == []:
            test_args = PAP_Optimal.pap_normal_data(input_normal)
            opt_test_normal = PAP_Optimal.optimal_test(test_args, .05)

            power_message = f"Expected power: {opt_test_normal['power'].round(2)}"

            t = opt_test_normal["t"].round(2)
            t_filtered = t.query('t > 0').sort_values('t')        
        else:
            opt_simple_normal = PAP_Simple.pap_normal_data_simple(input_normal, .05)
            power_message = f"Expected power: {opt_simple_normal['power'].round(2)}"
            t = opt_simple_normal["test"]     
            t_filtered = t
        
        table = dash_table.DataTable(
                data=t_filtered.to_dict('records'),
                columns=[{'name': i, 'id': i} for i in t_filtered.columns],
                style_cell={'minWidth': '100px', 'maxWidth': '200px'},
            )  
        
        return power_message, table, t.to_csv(index=False)

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