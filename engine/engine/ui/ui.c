#include "ui.h"

Status ui_init(UI* ui, Canvas* canvas)
{
	Status status = font_init(&ui->font);
	if (STATUS_OK != status)
	{
		return status;
	}

	ui->canvas = canvas;

	return STATUS_OK;
}

void ui_draw(UI* ui, float upscaling_factor)
{
	for (int i = 0; i < ui->text_count; ++i)
	{
		// TODO: Text needs to be scaled. Or how do I do it not scaled.
		text_draw(ui->canvas, &ui->text[i], &ui->font, upscaling_factor);
	}
}

void ui_destroy(UI* ui)
{

}
