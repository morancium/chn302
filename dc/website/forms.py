from django import forms
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit

class ParameterForm(forms.Form):

    feed_conc = forms.FloatField(
        label = "Feed Concentration",
        help_text = "Concentration of feed component",
    )
    feed_loc = forms.IntegerField(
        label = "Feed Location",
        help_text = "Tray where feed is inserted",
    )
    feed_flow = forms.FloatField(
        label = "Feed Flow Rate",
        help_text = "Feed flowrate",
    )
    feed_temp = forms.FloatField(
        label = "Feed Temperature",
        help_text = "Temperature of the feed",
    )
    pressure = forms.FloatField(
        label = "Pressure",
        help_text = "Pressure of the column",
    )
    distillate_flow = forms.FloatField(
        label = "Distillate Flow Rate",
        help_text = "Distilate flowrate",
    )
    tray_count = forms.IntegerField(
        label = "Number of Trays",
        help_text = "Number of trays",
    )
    reflux_ratio = forms.FloatField(
        label = "Reflux Ratio",
        help_text = "Reflux Ratio",
    )
    condenser_type = forms.ChoiceField(
        choices= ((0, "Partial"), (1, "Total"),),
        label = "Condenser Type",
        help_text = "Type of condenser",

    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_id = 'id-exampleForm'
        self.helper.form_class = 'blueForms'
        self.helper.form_method = 'post'
        self.helper.form_action = ''

        self.helper.add_input(Submit('submit', 'Submit'))