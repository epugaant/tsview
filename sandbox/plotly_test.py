import plotly.graph_objects as go

fig = go.Figure()

fig.add_trace(go.Scatter(x=[1, 2, 3, 4], y=[10, 15, 13, 17], 
                         error_y=dict(
                                type='data', # value of error bar given in data coordinates
                                array=[0.02, 0.03, 0.1, 0.2],
                                visible=True),
                         name='filter 1'))
fig.add_trace(go.Scatter(x=[1, 2, 3, 4], y=[16, 5, 11, 9],  
                         error_y=dict(
                                type='data', # value of error bar given in data coordinates
                                array=None,
                                visible=True),
                         name='filter 2'))

print(str(fig.to_json()))
fig.update_layout(legend_title_text = "Mission")
fig.update_xaxes(title_text="x")
fig.update_yaxes(title_text="y")
fig.show()
pass