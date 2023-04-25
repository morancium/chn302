from django.shortcuts import render, redirect
from django.conf import settings
from .forms import ParameterForm
from django.http import HttpResponseRedirect
from utils import model

# Create your views here.
def index(request):

    if request.method == "POST":
        form = ParameterForm(request.POST)
        if form.is_valid():
            request.session['data'] = form.cleaned_data
            return HttpResponseRedirect("results/")
    else:
        context ={}
        context['form']= ParameterForm()

        return render(request, str(settings.BASE_DIR) + "/static/templates/index.html", context)

def results(request):
    context = {}
    context['form_data'] = request.session['data']
    feed_conc = context['form_data']['feed_conc']
    feed_loc = context['form_data']['feed_loc']
    feed_flow = context['form_data']['feed_flow']
    feed_temp = context['form_data']['feed_temp']
    pressure = context['form_data']['pressure']
    distillate_flow = context['form_data']['distillate_flow']
    tray_count = context['form_data']['tray_count']
    reflux_ratio = context['form_data']['reflux_ratio']
    condenser_type = context['form_data']['condenser_type']
    
    return render(request, str(settings.BASE_DIR) + "/static/templates/results.html", context)